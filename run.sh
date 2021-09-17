#!/bin/bash

HDFS_OUTPUT_DIR="hdfs://bic07/user/ybw/ssd/output"
HDFS_DATASET_DIR="hdfs://bic07/user/ybw/ssd/Dataset"
RESULT_DIR=/home/ybw/Reproduce/thrill/result
THRILL_DIR=/home/ybw/Reproduce/thrill
THRILL_BUILD_DIR=${THRILL_DIR}/original-build
HOSTS_LIST=("bic01" "bic01 bic02" "bic01 bic02 bic03 bic04" "bic01 bic02 bic03 bic04 bic05 bic07 bic08 bic09")
HOST_NUM=(1 2 4 8)
CID=3
HN=${HOST_NUM[${CID}]}
HOSTS=${HOSTS_LIST[${CID}]}
echo ${HN}
echo ${HOSTS}

hdfs dfs -rm -r -skipTrash ${HDFS_OUTPUT_DIR}/thrill_pagerank
hdfs dfs -mkdir -p ${HDFS_OUTPUT_DIR}/thrill_pagerank
${THRILL_DIR}/run/ssh/invoke.sh -x LIBHDFS3_CONF=/etc/hadoop/conf/hdfs-site.xml -h "${HOSTS}" -w 28 ${THRILL_BUILD_DIR}/examples/page_rank/page_rank_run_end ${HDFS_DATASET_DIR}/twitter-2010.text -n 10 -j true -o ${HDFS_OUTPUT_DIR}/thrill_pagerank/result 1> ${RESULT_DIR}/pr_8_ssd.stdout.txt 2> ${RESULT_DIR}/pr_8_ssd.stderr.txt

${THRILL_DIR}/run/ssh/invoke.sh -x LIBHDFS3_CONF=/etc/hadoop/conf/hdfs-site.xml -h "${HOSTS}" -w 28 ${THRILL_BUILD_DIR}/examples/terasort/CC ${HDFS_DATASET_DIR}/twitter-2010.text 1> ${RESULT_DIR}/cc_8_ssd.stdout.txt 2> ${RESULT_DIR}/cc_8_ssd.stderr.txt

hdfs dfs -rm -r ssd/output/wordcount
hdfs dfs -mkdir -p ssd/output/wordcount
${THRILL_DIR}/run/ssh/invoke.sh -x LIBHDFS3_CONF=/etc/hadoop/conf/hdfs-site.xml -h "${HOSTS}" -w 28 ${THRILL_BUILD_DIR}/examples/word_count/word_count_run_end ${HDFS_DATASET_DIR}/deepmind-gutenberg.text -o ${HDFS_OUTPUT_DIR}/wordcount/result 1> ${RESULT_DIR}/wc_8_ssd.stdout.txt 2> ${RESULT_DIR}/wc_8_ssd.stderr.txt

hdfs dfs -rm -r ssd/output/terasort
hdfs dfs -mkdir ssd/output/terasort
${THRILL_DIR}/run/ssh/invoke.sh -x LIBHDFS3_CONF=/etc/hadoop/conf/hdfs-site.xml -h "${HOSTS}" -w 28 ${THRILL_BUILD_DIR}/examples/terasort/terasort ${HDFS_DATASET_DIR}/terasort-250G -s true -o ${HDFS_OUTPUT_DIR}/terasort/result 1> ${RESULT_DIR}/ts_8_ssd.stdout.txt 2> ${RESULT_DIR}/ts_8_ssd.stderr.txt

${THRILL_DIR}/run/ssh/invoke.sh -x LIBHDFS3_CONF=/etc/hadoop/conf/hdfs-site.xml -h "${HOSTS}" -w 28 ${THRILL_BUILD_DIR}/examples/k-means/k-means_run -n 10 64 500 hdfs://bic07/user/lixiaohan/data-kmeans/data-Kmeans-80000000-64 1> ${RESULT_DIR}/km_8_ssd.stdout.txt 2> ${RESULT_DIR}/km_8_ssd.stderr.txt

${THRILL_DIR}/run/ssh/invoke.sh -x THRILL_RAM=808527593472 -x LIBHDFS3_CONF=/etc/hadoop/conf/hdfs-site.xml -h "${HOSTS}" -w 28 ${THRILL_BUILD_DIR}/examples/logistic_regression/logistic_regression_run ${HDFS_DATASET_DIR}/MNIST01-12665-10k-784.text -n 10 -g 1.0 -e 0.0000001 1> ${RESULT_DIR}/lr_8_ssd.stdout.txt 2> ${RESULT_DIR}/lr_8_ssd.stderr.txt

