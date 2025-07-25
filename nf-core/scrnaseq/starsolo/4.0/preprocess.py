from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Make a wide samplesheet with the columns
# sample, fastq_1, fastq_1
samplesheet = (
    ds.files
    .assign(
        readType=lambda d: d.apply(
            lambda r: r.get("readType", "R"),
            axis=1
        ),
        lane=lambda d: d.apply(
            lambda r: r.get("lane", 1),
            axis=1
        ),
    )
    .query("readType == 'R'")
    .pivot(
        index=["sampleIndex", "lane", "sample", "dataset"],
        columns="read",
        values="file"
    )
    .rename(columns=lambda i: f"fastq_{int(i)}")
    .reset_index()
    .reindex(columns=["sample", "fastq_1", "fastq_2"])
)

ds.logger.info("Formatted samplesheet:")
ds.logger.info(samplesheet.to_csv(index=None))
assert samplesheet.shape[0] > 0, "No FASTQ files detected -- there may be an error with user selection or data ingest"

# Write out to a file
samplesheet.to_csv("samplesheet.csv", index=None)

# Flatten the file selection params
for kw in ["star_index", "gtf"]:
    if kw in ds.params and len(ds.params[kw]) > 0:
        ds.add_param(
            kw,
            ds.params[kw][0],
            overwrite=True
        )
        ds.logger.info(f"Formatted parameter {kw}: {ds.params[kw]}")
    else:
        raise ValueError(f"Missing required parameter: {kw} in alignment options.")

# Remove the "/SA" from the STAR index path
assert ds.params["star_index"].endswith("/SA"), "Expected STAR index path to end with '/SA'"
ds.add_param(
    "star_index",
    ds.params["star_index"].replace("/SA", ""),
    overwrite=True
)
ds.logger.info(f"Updated STAR index path: {ds.params['star_index']}")
