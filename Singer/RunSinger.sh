
vcf=$1
ratio=$2
mkdir ${vcf}_${ratio}
path_to_singer=/home/jobste/github/SINGER/singer-0.1.7-beta-linux-x86_64/
${path_to_singer}/parallel_singer -Ne 6.22e5 -m 3.4e-9 \
-vcf ${vcf} -output ${vcf}_${ratio}/${vcf}_Singer \
-n 100 -thin 10 -polar 0.5 -ratio ${ratio} -num_cores 4
