import allel
from multiprocessing import Process

vcfpath='/home/dennist/lstm_data/cease/variants_bycohort/combined_cohorts/'
vcfs = [f'{vcfpath}CM023248.noindels.filtered.ann.vcf.gz',f'{vcfpath}CM023249.noindels.filtered.ann.vcf.gz',f'{vcfpath}CM023249.noindels.filtered.ann.vcf.gz']
zarrpath = '/home/dennist/lstm_data/cease/variants_bycohort/combined_cohorts/annotated_zarr'

def chr2():
    print('converting chr2')
    allel.vcf_to_zarr('/home/dennist/lstm_data/cease/variants_bycohort/combined_cohorts/vcf/CM023248.noindels.filtered.ann.vcf.gz', f'{zarrpath}/annotations.zarr', fields=['variants/ANN', 'variants/POS', 'calldata/GT'], group='CM023248', transformers=allel.ANNTransformer(), exclude_fields = ['variants/LOF', 'variants/NMD'], overwrite=True)

def chr3():
    print('converting chr3')
    allel.vcf_to_zarr('/home/dennist/lstm_data/cease/variants_bycohort/combined_cohorts/vcf/CM023249.noindels.filtered.ann.vcf.gz', f'{zarrpath}/annotations.zarr', fields=['variants/ANN', 'variants/POS', 'calldata/GT'], group='CM023249', transformers=allel.ANNTransformer(), exclude_fields = ['variants/LOF', 'variants/NMD'], overwrite=True)

def chrX():
    print('converting chrX')
    allel.vcf_to_zarr('/home/dennist/lstm_data/cease/variants_bycohort/combined_cohorts/vcf/CM023250.noindels.filtered.ann.vcf.gz', f'{zarrpath}/annotations.zarr', fields=['variants/ANN', 'variants/POS', 'calldata/GT'], group='CM023250', transformers=allel.ANNTransformer(), exclude_fields = ['variants/LOF', 'variants/NMD'], overwrite=True)

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
