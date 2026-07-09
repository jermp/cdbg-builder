GGCAT
===

	/usr/bin/time -v ggcat build -k 31 -j 32 -l /mnt/hd2/pibiri/DNA/filenames/blackwell-100k_filenames.txt -s 1 -c -o /mnt/hd2/pibiri/DNA/bw-100k-ggcat -m 16 -t /mnt/hd2/pibiri/DNA/tmp_dir

	/usr/bin/time -v ggcat build -k 31 -j 48 -l /mnt/hd2/pibiri/DNA/filenames/blackwell-661k_filenames.txt -s 1 -c -o /mnt/hd2/pibiri/DNA/bw-661k-ggcat -m 64 -t /mnt/hd2/pibiri/DNA/tmp_dir

	/usr/bin/time -v ggcat build -k 31 -j 32 -l /mnt/hd2/pibiri/DNA/filenames/se_4546_filenames.txt -s 1 -c -o /mnt/hd2/pibiri/DNA/se-4546-ggcat -m 8 -t /mnt/hd2/pibiri/DNA/tmp_dir

CF3
===

	export PARLAY_NUM_THREADS=32
	ulimit -n 32768
	source ~/my_local_env/bin/activate
	/usr/bin/time -v ./cuttlefish build -l /mnt/hd2/pibiri/DNA/filenames/blackwell-100k_filenames.txt -k 31 -o /mnt/hd2/pibiri/DNA/cf-bw-100k -w /mnt/hd2/pibiri/DNA/tmp_dir --ref --color -c 1