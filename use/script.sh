#!/bin/bash
index_array=("new_request_cnt_throughputs" \
             "num_of_node_throughputs" "new_request_cnt_use_memory_ratio" "new_request_cnt_use_channel_ratio" \
             "new_request_cnt_runtime"  \)
#index_array=("min_fidelity_encode_ratio") 


for i in ${index_array[@]}
do
    #rm $i/*eps
    python3 $i/ChartGenerator.py
    cp $i/*.eps eps/QEC_$i.eps 
    cp $i/*.jpg jpg/QEC_$i.jpg 
    
done