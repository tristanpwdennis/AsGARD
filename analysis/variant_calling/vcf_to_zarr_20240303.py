import allel
from multiprocessing import Process
from sgkit.io.vcf import vcf_to_zarr
from dask.distributed import Client

vcfs = ['~/lstm_data/cease/variants_bycohort/combined_cohorts/vcf/CM023248.noindels.filtered.ann.vcf.gz','~/lstm_data/cease/variants_bycohort/combined_cohorts/vcf/CM023249.noindels.filtered.ann.vcf.gz','~/lstm_data/cease/variants_bycohort/combined_cohorts/vcf/CM023250.noindels.filtered.ann.vcf.gz']


def chr2():
	print('converting chr2')
	if __name__ == "__main__":

		client = Client(n_workers=10, threads_per_worker=1)

		zarr  = 'zarr/combined_cohorts.CM023248.zarr'

		vcf_to_zarr(vcfs[0],
			zarr,
			max_alt_alleles=1,
			#exclude_fields="INFO/ANN",
            fields=["FORMAT/GT", "FORMAT/DP", "FORMAT/GQ", "FORMAT/AD"],
            tempdir='zarr/tmp/')


def chr3():
	print('converting chr3')
	if __name__ == "__main__":

		client = Client(n_workers=10, threads_per_worker=1)

		zarr  = 'zarr/combined_cohorts.CM023249.zarr'

		vcf_to_zarr(vcfs[1],
			zarr,
			max_alt_alleles=1,
			#exclude_fields="INFO/ANN",
            fields=["FORMAT/GT", "FORMAT/DP", "FORMAT/GQ", "FORMAT/AD"],
            tempdir='zarr/tmp/')

def chrX():
	print('converting chrX')
	if __name__ == "__main__":

		client = Client(n_workers=10, threads_per_worker=1)

		zarr  = 'zarr/combined_cohorts.CM023250.zarr'

		vcf_to_zarr(vcfs[2],
			zarr,
			max_alt_alleles=1,
			#exclude_fields="INFO/ANN",
            fields=["FORMAT/GT", "FORMAT/DP", "FORMAT/GQ", "FORMAT/AD"],
            tempdir='zarr/tmp/')



#do in parallel
if __name__ == "__main__":
    p1 = Process(target=chr2)
    p1.start()
    p2 = Process(target=chr3)
    p2.start()
    p3 = Process(target=chrX)
    p3.start()

    p1.join()
    p2.join()
    p3.join()
