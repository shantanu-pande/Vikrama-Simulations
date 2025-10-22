from co2_ess_sim import run_sim, plot_results


def main() -> None:
	# Run a 24-hour simulation for a 100 MW block and plot results
	df, summary = run_sim(hours=24, dt_s=60, target_power_MW=100.0, fluid="CO2")
	print("Technical feasibility snapshot")
	print(f"- Round-trip efficiency (sim): {summary['RTE']*100:.1f}%")
	print(f"- Estimated block capacity (model): ~{summary['capacity_est_MWh']:.1f} MWh")
	print(f"- Illustrative LCoS (prelim.): Rs {summary['LCOS_Rs_per_kWh']:.2f} / kWh")
	plot_results(df, "CO₂ Compression ESS — 24h Operation")


if __name__ == "__main__":
	main()
