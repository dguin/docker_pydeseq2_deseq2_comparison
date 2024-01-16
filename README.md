<h1>Framework for comparing pyDESeq2 vs. DESeq2</h1>

This is a general framework for quickly running a comparison between pyDESeq2 and DESeq2 stable. The comparison can be run fully inside the provided docker container without installing any dependencies locally. To run the comparison in the docker follow the following steps:
First, build the docker image. 
Note: The R docker containers may perform poorly in the new MAC ARM architecture. If so build and use the docker using the `--platform linux/amd64` flag.

To build the docker from the root do:
```sh
docker build \
  --build-arg PYDESEQ2_VERSION={PYDESEQ2_VERSION} \
  -t {DOCKER_IMAGE_NAME} .
```
For example, to compare pyDESeq2 0.4.4 to R DESeq2 do:
```sh
docker build \
  --build-arg PYDESEQ2_VERSION=0.4.4 \
  -t {DOCKER_IMAGE_NAME} .
```
with your choice of the docker image name.

Before running the docker container build the count and design files and save them to a director of your choice as csvs. The count file should have genes as the index and sample IDs as column names in the csv. The design csv file must contain the sample IDs or column names from the count csv as the index. In addition it must contain a condition column identifying either `control` or `test` for each sample ID. Additional columns can be provided containing any other variables to identify any secondary effects. For the comparison, DE results will be returned by comparing the test samples to the control. Save the count file as `counts.csv` and the design file as `metadata.csv` to the same directory. For an example of both `counts.csv` and metadata.csv see example_data.

Once counts have been saved the docker container can be run using
```sh
docker run -it \
  -v ${DATA_PATH}:/opt/data \
  -v ${OUTPUT_PATH}:/opt/output \
  {DOCKER_IMAGE_NAME}
```
**Note 1:** The path defined by `${DATA_PATH}` in the docker run command should be the same path where the `counts.csv` and `metadata.csv` have been saved.
**Note 2:** Both `${DATA_PATH}` and `${OUTPUT_PATH}` provided to the docker run -it command should be absolute paths.

The path `${OUTPUT_PATH}` contains the results of the comparison.

Once done the outputs will be available at the output path.
Output file list:
1. `{r,py}_results_post_shrunk.csv`: Fold change csv post-shrinkage.
2. `{r,py}_results_pre_shrunk.csv`: Fold change csv pre-shrinkage.
3. `{r,py}_size_factors.csv`: Size factors.
Results from R are prefixed with `r` and results from python are prefixed with `py`.