conda create -n accelerometer python=3.9 openjdk pip -y
conda activate accelerometer

pip install accelerometer
javac -cp java/JTransforms-3.1-with-dependencies.jar java/*.java


mkdir -p results

# Define the input file
RAW_FILE="/mnt/project/Bulk/Activity/Raw/10/1000056_90001_0_0.cwa"

# --- FIX 3: Copy the raw file to the local writable directory before processing ---
echo "Copying input file to local directory..."
cp "$RAW_FILE" .
LOCAL_RAW_FILE=$(basename "$RAW_FILE")

# Run the processing command on the local copy of the file
echo "--- Processing file: $LOCAL_RAW_FILE ---"


accProcess $RAW_FILE \
  --activityClassification False \
  --deleteIntermediateFiles False \
  --outputFolder results

for f in results/*-epoch.csv.gz; do
  id=${f%%-*}
  zcat "$f" \
  | awk -F',' -v id="$id" '
      NR==1 {
        for (i=1;i<=NF;i++){
          if ($i=="time") t=i; if ($i=="temp") p=i; if ($i=="enmoTrunc") e=i;
          if ($i=="dataErrors") d=i; if ($i=="samples") s=i;
          if ($i=="clipsBeforeCalibr") c1=i; if ($i=="clipsAfterCalibr") c2=i }
        print "id,time,temp,enmoTrunc,dataErrors,samples,clipsBefore,clipsAfter"; next
      }
      {print id","$t","$p","$e","$d","$s","$c1","$c2}
    '
done > temperature_epochs.csv
gzip temperature_epochs.csv

# 1.  Write the header row once
echo "id,goodWear,goodCal,calErrAfter_mg,wearDays,\
xOffset_g,yOffset_g,zOffset_g,\
clipsBefore,clipsAfter,startTime,endTime" \
> summary_covariates.csv

# 2.  Append one CSV line per summary.json
for jf in results/*-summary.json ; do
  id=$(basename "$jf" | cut -d'-' -f1)        # 1000056_90001_0_0
  jq -r --arg id "$id" '
    [ $id,
      ."quality-goodWearTime",
      ."quality-goodCalibration",
      ."calibration-errsAfter(mg)",
      ."wearTime-overall(days)",
      ."calibration-xOffset(g)",
      ."calibration-yOffset(g)",
      ."calibration-zOffset(g)",
      ."clipsBeforeCalibration",
      ."clipsAfterCalibration",
      ."file-startTime",
      ."file-endTime"
    ] | @csv
  ' "$jf" >> summary_covariates.csv
done

dx mkdir -p /temp_ukb_results

dx upload temperature_epochs.csv.gz \
  --path /temp_ukb_results/ \
  --wait            # optional: block until the file finishes closing

dx upload summary_covariates.csv \
  --path /temp_ukb_results/ \
  --wait