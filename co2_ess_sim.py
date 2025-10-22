"""
CO2/Air Compression-Based ESS — Simplified System Simulation

What this does
- Simulates charge/discharge of a CO₂/air compression energy storage plant
- Tracks tank pressure/temperature, TES temperature, SoC, power, energy
- Produces plots for quick visualization

Notes
- Uses CoolProp for real-gas CO₂ if available; otherwise ideal-gas fallback.
- Models are lumped/first-order and meant for concepting, not detailed design.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Try real-gas properties for CO2; fall back to ideal gas if missing
USE_COOLPROP = False
try:
    from CoolProp.CoolProp import PropsSI  # type: ignore
    USE_COOLPROP = True
except Exception:
    warnings.warn(
        "CoolProp not found — using ideal-gas approximations. Install with: pip install CoolProp"
    )


# -----------------------------
# Physical helpers
# -----------------------------
R_UNIVERSAL = 8.314462618  # J/mol-K
MOLAR_MASS = {
    "CO2": 44.0095e-3,  # kg/mol
    "Air": 28.97e-3,    # kg/mol
}


def cp_k_for(fluid: str, T: float) -> Tuple[float, float]:
    """Return rough cp and k for ideal-gas fallback (J/kg-K, ratio of specific heats)."""
    if fluid == "CO2":
        cp = 0.85e3  # in the 250–350 K band (approx.)
        k = 1.30
    else:  # Air
        cp = 1.005e3
        k = 1.40
    return cp, k


def eos_density(fluid: str, p: float, T: float) -> float:
    """Compute density (kg/m3)."""
    if USE_COOLPROP and fluid == "CO2":
        return float(PropsSI("D", "P", p, "T", T, "CO2"))
    # Ideal gas fallback
    R = R_UNIVERSAL / MOLAR_MASS.get(fluid, MOLAR_MASS["Air"])
    return p / (R * T)


def enthalpy(fluid: str, p: float, T: float) -> float:
    """Return specific enthalpy (J/kg). Only differences matter in this model."""
    if USE_COOLPROP and fluid == "CO2":
        return float(PropsSI("H", "P", p, "T", T, "CO2"))
    cp, _ = cp_k_for(fluid, T)
    return cp * (T - 273.15)  # referenced to 0C; only differences used


def entropy(fluid: str, p: float, T: float) -> float:
    """Return specific entropy (J/kg-K). Only differences used for isentropic calcs."""
    if USE_COOLPROP and fluid == "CO2":
        return float(PropsSI("S", "P", p, "T", T, "CO2"))
    cp, k = cp_k_for(fluid, T)
    R = R_UNIVERSAL / MOLAR_MASS.get(fluid, MOLAR_MASS["Air"])
    return cp * math.log(T / 298.15) - R * math.log(p / 1e5)


# -----------------------------
# Component models
# -----------------------------


@dataclass
class Compressor:
    eta_isentropic: float = 0.82
    stages: int = 3
    max_power_W: float = 100e6  # 100 MW limit

    def compress(self, fluid: str, m_dot: float, p_in: float, T_in: float, p_out: float) -> Tuple[float, float, float]:
        """Multi‑stage compression with intercooling to ambient between stages.

        Returns: (p_out, T_after_intercool, electrical_power_W)
        """
        if m_dot <= 0.0 or p_out <= p_in:
            return p_in, T_in, 0.0

        p_ratio_total = p_out / p_in
        p_ratio_stage = p_ratio_total ** (1.0 / self.stages)
        T_stage_in = T_in
        p_stage_in = p_in
        work_total = 0.0
        T_ambient_cooler = 298.15  # 25C intercooling target

        for _ in range(self.stages):
            p_stage_out = p_stage_in * p_ratio_stage
            w, T_stage_out = self._stage_work(fluid, p_stage_in, T_stage_in, p_stage_out)
            work_total += w * m_dot
            # Intercool to ambient (reduces next stage work)
            T_stage_in = T_ambient_cooler
            p_stage_in = p_stage_out

        # Cap by available power
        if work_total > self.max_power_W:
            scale = self.max_power_W / max(work_total, 1e-9)
            m_dot *= scale
            work_total = self.max_power_W

        return p_out, T_stage_in, work_total

    def _stage_work(self, fluid: str, p1: float, T1: float, p2: float) -> Tuple[float, float]:
        """Return (specific work J/kg, T_out K) for one stage."""
        if USE_COOLPROP and fluid == "CO2":
            s1 = entropy(fluid, p1, T1)
            T2s = float(PropsSI("T", "P", p2, "S", s1, "CO2"))
            h1 = enthalpy(fluid, p1, T1)
            h2s = enthalpy(fluid, p2, T2s)
            w_is = h2s - h1
            w = w_is / max(self.eta_isentropic, 1e-6)
            cp, _ = cp_k_for(fluid, (T1 + T2s) / 2)
            T2 = T1 + w / cp
            return w, T2
        else:
            cp, k = cp_k_for(fluid, T1)
            T2s = T1 * (p2 / p1) ** ((k - 1) / k)
            w_is = cp * (T2s - T1)
            w = w_is / max(self.eta_isentropic, 1e-6)
            T2 = T1 + w / cp
            return w, T2


@dataclass
class Turbine:
    eta_isentropic: float = 0.86
    max_power_W: float = 100e6  # 100 MW

    def expand(
        self,
        fluid: str,
        m_dot: float,
        p_in: float,
        T_in: float,
        p_out: float,
        T_reheat: Optional[float] = None,
    ) -> Tuple[float, float, float]:
        """Expand through turbine with optional reheating upstream of turbine."""
        if m_dot <= 0.0 or p_in <= p_out:
            return p_in, T_in, 0.0

        T1 = T_in if T_reheat is None else max(T_in, T_reheat)
        if USE_COOLPROP and fluid == "CO2":
            s1 = entropy(fluid, p_in, T1)
            T2s = float(PropsSI("T", "P", p_out, "S", s1, "CO2"))
            h1 = enthalpy(fluid, p_in, T1)
            h2s = enthalpy(fluid, p_out, T2s)
            w_is = h1 - h2s
            w = self.eta_isentropic * w_is
            power = w * m_dot
            if power > self.max_power_W:
                scale = self.max_power_W / max(power, 1e-9)
                m_dot *= scale
                power = self.max_power_W
            cp, _ = cp_k_for(fluid, (T1 + T2s) / 2)
            T2 = T1 - w / cp
            return p_out, T2, power
        else:
            cp, k = cp_k_for(fluid, T1)
            T2s = T1 * (p_out / p_in) ** ((k - 1) / k)
            w_is = cp * (T1 - T2s)
            w = self.eta_isentropic * w_is
            power = w * m_dot
            if power > self.max_power_W:
                scale = self.max_power_W / max(power, 1e-9)
                m_dot *= scale
                power = self.max_power_W
            T2 = T1 - w / cp
            return p_out, T2, power


@dataclass
class ThermalStorage:
    # Lumped TES: single-node temperature with heat loss to ambient
    mass_kg: float = 5_000_000.0
    heat_capacity_J_per_kgK: float = 900.0
    UA_W_per_K: float = 50_000.0
    T_init_K: float = 320.0
    T_min_K: float = 290.0
    T_max_K: float = 650.0

    def __post_init__(self) -> None:
        self.T = self.T_init_K

    @property
    def C(self) -> float:
        return self.mass_kg * self.heat_capacity_J_per_kgK

    def add_heat(self, Q_J: float, dt_s: float, T_amb_K: float) -> None:
        loss = self.UA_W_per_K * (self.T - T_amb_K) * dt_s
        dE = Q_J - loss
        self.T = float(np.clip(self.T + dE / max(self.C, 1e-9), self.T_min_K, self.T_max_K))

    def withdraw_heat(self, Q_req_J: float, dt_s: float, T_amb_K: float) -> float:
        loss = self.UA_W_per_K * (self.T - T_amb_K) * dt_s
        available = max((self.T - self.T_min_K) * self.C - loss, 0.0)
        Q_out = min(Q_req_J, available)
        dE = -Q_out - loss
        self.T = float(np.clip(self.T + dE / max(self.C, 1e-9), self.T_min_K, self.T_max_K))
        return Q_out


@dataclass
class StorageTank:
    volume_m3: float = 2000.0
    p_min_Pa: float = 70e5
    p_max_Pa: float = 90e5
    T_init_K: float = 273.15 + (-5)
    UA_W_per_K: float = 30_000.0
    fluid: str = "CO2"

    def __post_init__(self) -> None:
        self.p = 0.5 * (self.p_min_Pa + self.p_max_Pa)
        self.T = self.T_init_K
        self.m = eos_density(self.fluid, self.p, self.T) * self.volume_m3

    def update_state_from_m(self) -> None:
        if USE_COOLPROP and self.fluid == "CO2":
            rho = self.m / self.volume_m3
            self.p = float(PropsSI("P", "D", rho, "T", self.T, "CO2"))
        else:
            R = R_UNIVERSAL / MOLAR_MASS.get(self.fluid, MOLAR_MASS["Air"])
            self.p = self.m / self.volume_m3 * R * self.T

    def add_mass(self, m_in: float, T_in_K: float, dt_s: float, T_amb_K: float) -> None:
        if m_in > 0.0:
            self.T = (self.m * self.T + m_in * T_in_K) / (self.m + m_in + 1e-9)
            self.m += m_in
        loss = self.UA_W_per_K * (self.T - T_amb_K) * dt_s
        cp, _ = cp_k_for(self.fluid, self.T)
        self.T -= loss / max(cp * self.m, 1e-9)
        self.update_state_from_m()

    def remove_mass(self, m_out: float, dt_s: float, T_amb_K: float) -> None:
        self.m = max(self.m - m_out, 1e-6)
        loss = self.UA_W_per_K * (self.T - T_amb_K) * dt_s
        cp, _ = cp_k_for(self.fluid, self.T)
        self.T -= loss / max(cp * self.m, 1e-9)
        self.update_state_from_m()

    @property
    def soc(self) -> float:
        return float(np.clip((self.p - self.p_min_Pa) / (self.p_max_Pa - self.p_min_Pa), 0.0, 1.0))


@dataclass
class Plant:
    fluid: str = "CO2"
    compressor: Compressor = Compressor()
    turbine: Turbine = Turbine()
    tank: StorageTank = StorageTank()
    tes: ThermalStorage = ThermalStorage()
    T_ambient_K: float = 298.15

    def step(self, mode: str, P_target_W: float, dt_s: float) -> Dict[str, float]:
        """Advance plant by one time step.

        mode: 'charge' | 'discharge' | 'idle'
        Returns dict of key signals for logging/plots.
        """
        m_dot = 0.0
        P_elec = 0.0
        Q_to_TES = 0.0
        Q_from_TES = 0.0

        if mode == "charge" and self.tank.p < self.tank.p_max_Pa * 0.999:
            m_dot = 5.0  # initial guess kg/s
            for _ in range(10):
                p_out = min(self.tank.p_max_Pa, self.tank.p * 1.1)
                p_out = max(p_out, self.tank.p + 1e4)
                _, T_cool_to_tank, P_in = self.compressor.compress(self.fluid, m_dot, self.tank.p, self.T_ambient_K, p_out)
                if P_in == 0:
                    break
                scale = float(np.clip(P_target_W / (P_in + 1e-9), 0.2, 5.0))
                m_dot *= scale
            # Final apply
            _, T_cool_to_tank, P_in = self.compressor.compress(self.fluid, m_dot, self.tank.p, self.T_ambient_K, p_out)
            P_elec = P_in
            Q_to_TES = 0.9 * P_in * dt_s  # send ~90% of compressor work to TES
            self.tes.add_heat(Q_to_TES, dt_s, self.T_ambient_K)
            self.tank.add_mass(m_dot * dt_s, T_cool_to_tank, dt_s, self.T_ambient_K)

        elif mode == "discharge" and self.tank.p > self.tank.p_min_Pa * 1.001:
            T_reheat = min(self.tes.T, 600.0)
            m_dot = 5.0
            for _ in range(10):
                _, _, P_out = self.turbine.expand(self.fluid, m_dot, self.tank.p, self.tank.T, max(self.tank.p_min_Pa, self.tank.p * 0.9), T_reheat=T_reheat)
                if P_out == 0:
                    break
                scale = float(np.clip(P_target_W / (P_out + 1e-9), 0.2, 5.0))
                m_dot *= scale
            _, T_out, P_out = self.turbine.expand(self.fluid, m_dot, self.tank.p, self.tank.T, max(self.tank.p_min_Pa, self.tank.p * 0.9), T_reheat=T_reheat)
            P_elec = -P_out  # negative means generation
            if T_reheat > self.tank.T:
                cp, _ = cp_k_for(self.fluid, (T_reheat + self.tank.T) / 2)
                Q_req = cp * m_dot * (T_reheat - self.tank.T) * dt_s
                Q_from_TES = self.tes.withdraw_heat(Q_req, dt_s, self.T_ambient_K)
            self.tank.remove_mass(m_dot * dt_s, dt_s, self.T_ambient_K)

        else:
            # Idle losses only
            self.tank.add_mass(0.0, self.tank.T, dt_s, self.T_ambient_K)
            self.tes.add_heat(0.0, dt_s, self.T_ambient_K)

        return dict(
            p_Pa=self.tank.p,
            T_tank_K=self.tank.T,
            T_tes_K=self.tes.T,
            soc=self.tank.soc,
            P_elec_W=P_elec,
            Q_to_TES_J=Q_to_TES,
            Q_from_TES_J=Q_from_TES,
            m_dot_kg_s=m_dot,
        )


def synthetic_profiles(hours: int = 24, dt_s: int = 60) -> pd.DataFrame:
    """Create simple renewable and demand profiles to drive EMS decisions."""
    t = np.arange(0, hours * 3600 + dt_s, dt_s)
    # Renewable profile: solar-like + wind variation (per unit)
    solar = np.clip(np.sin((t / 3600 - 6) * np.pi / 12), 0, 1)
    wind = 0.4 + 0.2 * np.sin(t / 3600 * 2 * np.pi / 6 + 1.0)
    re_pu = np.clip(0.6 * solar + 0.4 * wind, 0, 1.2)

    # Demand profile: morning/evening peaks (per unit)
    demand = 0.7 + 0.2 * np.sin((t / 3600 - 7) * np.pi / 12) + 0.2 * np.sin((t / 3600 - 19) * np.pi / 8)
    demand = np.clip(demand, 0.5, 1.2)

    rating_MW = 100.0
    re_MW = re_pu * rating_MW
    demand_MW = demand * rating_MW
    net_MW = re_MW - demand_MW

    df = pd.DataFrame(dict(t_s=t, re_MW=re_MW, demand_MW=demand_MW, net_MW=net_MW))
    return df


def run_sim(hours: int = 24, dt_s: int = 60, target_power_MW: float = 100.0, fluid: str = "CO2") -> Tuple[pd.DataFrame, Dict[str, float]]:
    """Run the simulation for given horizon and return (timeseries, summary)."""
    plant = Plant(fluid=fluid)
    prof = synthetic_profiles(hours=hours, dt_s=dt_s)
    out: List[Dict[str, float]] = []
    for _, row in prof.iterrows():
        net = float(row["net_MW"])
        if net > 5.0 and plant.tank.soc < 0.99:
            mode = "charge"
            P_target = min(target_power_MW, net) * 1e6
        elif net < -5.0 and plant.tank.soc > 0.01:
            mode = "discharge"
            P_target = min(target_power_MW, -net) * 1e6
        else:
            mode = "idle"
            P_target = 0.0
        res = plant.step(mode, P_target, dt_s)
        res.update(dict(
            t_s=float(row["t_s"]),
            mode=mode,
            re_MW=float(row["re_MW"]),
            demand_MW=float(row["demand_MW"]),
            net_MW=float(row["net_MW"]),
        ))
        out.append(res)

    df = pd.DataFrame(out)

    # Energy accounting
    df["P_MW"] = df["P_elec_W"] / 1e6
    df["E_in_MWh"] = np.where(df["P_MW"] > 0, df["P_MW"] * dt_s / 3600.0, 0.0)
    df["E_out_MWh"] = np.where(df["P_MW"] < 0, -df["P_MW"] * dt_s / 3600.0, 0.0)
    E_in = float(df["E_in_MWh"].sum())
    E_out = float(df["E_out_MWh"].sum())
    rte = (E_out / E_in) if E_in > 0 else 0.0

    # Simple LCoS sketch (illustrative only)
    capex_per_kWh = 8000.0  # Rs/kWh placeholder
    opex_rs_per_kWyr = 800.0
    capacity_MWh = 0.12 * plant.tank.volume_m3  # indicative
    power_MW = 100.0
    years = 20
    total_delivered_MWh = E_out if E_out > 0 else 1.0
    ann_capex_rs = capacity_MWh * 1e3 * capex_per_kWh / years
    ann_opex_rs = power_MW * 1e3 * opex_rs_per_kWyr
    # Approximate average energy per day scaled to year (very rough):
    daily_delivered_MWh = E_out
    lcos_rs_per_kWh = (ann_capex_rs + ann_opex_rs) / max(daily_delivered_MWh * 365, 1e-6)

    summary = dict(
        E_in_MWh=E_in,
        E_out_MWh=E_out,
        RTE=rte,
        LCOS_Rs_per_kWh=lcos_rs_per_kWh,
        capacity_est_MWh=capacity_MWh,
    )
    return df, summary


def plot_results(df: pd.DataFrame, title: str = "CO2/ Air Compression ESS Simulation") -> None:
    t_h = df["t_s"] / 3600.0
    fig, axs = plt.subplots(4, 1, figsize=(12, 14), sharex=True)

    axs[0].plot(t_h, df["re_MW"], label="Renewables (MW)", alpha=0.8)
    axs[0].plot(t_h, df["demand_MW"], label="Demand (MW)", alpha=0.8)
    axs[0].plot(t_h, df["P_MW"], label="Plant Power (+charge / -discharge)", lw=2)
    axs[0].legend()
    axs[0].set_ylabel("Power (MW)")
    axs[0].grid(True, alpha=0.3)

    axs[1].plot(t_h, df["soc"] * 100, label="State of Charge (%)", color="tab:green")
    axs[1].set_ylabel("SoC (%)")
    axs[1].grid(True, alpha=0.3)

    axs[2].plot(t_h, df["p_Pa"] / 1e5, label="Tank Pressure (bar)", color="tab:red")
    axs[2].plot(t_h, df["T_tank_K"] - 273.15, label="Tank Temp (C)", color="tab:orange", alpha=0.7)
    axs[2].legend()
    axs[2].set_ylabel("Pressure (bar) / Temp (C)")
    axs[2].grid(True, alpha=0.3)

    axs[3].plot(t_h, df["T_tes_K"] - 273.15, label="TES Temp (C)", color="tab:purple")
    axs[3].plot(t_h, df["Q_to_TES_J"].cumsum() / 3.6e9, label="Heat -> TES (MWh_th)")
    axs[3].plot(t_h, df["Q_from_TES_J"].cumsum() / 3.6e9, label="Heat <- TES (MWh_th)")
    axs[3].legend()
    axs[3].set_xlabel("Time (hours)")
    axs[3].set_ylabel("Temp (C) / Thermal Energy (MWh_th)")
    axs[3].grid(True, alpha=0.3)

    fig.suptitle(title)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    df, summary = run_sim(hours=24, dt_s=60, target_power_MW=100.0, fluid="CO2")
    print("Summary:", summary)
    plot_results(df, "CO₂ Compression ESS — 24h Operation")
