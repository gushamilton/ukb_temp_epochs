Data-dictionary
(temperature-epoch Parquet + summary-covariate CSV)

Column	Type	Units	Definition & QC notes
Temperature-epoch Parquet temperature_epochs_<RUN_NAME>.parquet (Zstd-compressed)			
id	string	–	Device/participant identifier taken from the raw filename (e.g. 1000056_90001_0_0). Primary join key across all tables.
time	timestamp (ns, tz = Europe/London)	–	Start of the 30-s epoch, exactly as emitted by accProcess. Use for day-night aggregation or cosinor fitting.
temp	double	°C	Mean calibrated wrist-sensor temperature over the 30 s window. Calibration uses the axis slopes/offsets stored in the summary JSON.
enmoTrunc	double	milli-g (mg)	Activity intensity (Euclidean Norm Minus One, negatives truncated to 0). Include as an acute-activity covariate.
dataErrors	int32	count	Raw 100 Hz packets missing/corrupt inside the epoch. QC: drop or down-weight rows with dataErrors > 0 (rare).
samples	int32	count	100 Hz samples in the epoch (normally 3000). Lets you weight means if partial epochs are kept.
clipsBefore	int32	count	ADC saturations (±8 g) before calibration. High values (>500) may indicate damaged sensor housing.
clipsAfter	int32	count	ADC saturations after calibration (should be ≤ clipsBefore). Persistent non-zero → suspect device fault.
Summary-covariate CSV summary_covariates_<RUN_NAME>.csv.gz			
id	string	–	Same identifier as above.
goodWear	0 / 1	–	≥ 72 h good wear (quality-goodWearTime). Recommended filter: keep 1.
goodCal	0 / 1	–	Autocalibration converged & residual error ≤ 10 mg. Filter: keep 1.
calErrAfter_mg	float	mg	RMS calibration residual (calibration-errsAfter(mg)). QC: exclude if > 10 mg or add as covariate.
wearDays	float	days	Total valid wear days (wearTime-overall(days)). Use as covariate to control for data density.
xOffset_g, yOffset_g, zOffset_g	float	g	Axis-specific bias terms from sphere-fit calibration. Optional device-bias covariates; drop if they don’t change results.
clipsBefore, clipsAfter	int	count	File-level sums of saturated readings. High counts suggest impacts or cracked casing; often filtered at >500.
startTime, endTime	ISO 8601 string	–	First & last sample timestamps. Useful to merge with season, ambient weather, or assessment-centre date.

Recommended QC workflow
Sample-level filters

keep rows where goodWear == 1, goodCal == 1, and calErrAfter_mg ≤ 10.

optionally drop devices with clipsAfter > 500 or temperature range outside 25 – 45 °C.

Epoch-level filters

discard epochs where dataErrors > 0.

weight by samples if you keep partial windows.