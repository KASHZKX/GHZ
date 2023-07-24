#!/bin/bash
index_array=("min_fidelity_throughputs" "min_fidelity_fail-finished_ratio" "new_request_cnt_throughputs" \
             "num_of_node_throughputs" "new_request_cnt_use_memory_ratio" "new_request_cnt_use_channel_ratio" \
             "new_request_cnt_runtime" "min_fidelity_encode_rate" "path_length_encode_cnt" "encode_path_ratio" \
             "new_request_cnt_waiting_time" "min_fidelity_waiting_time" "encode_math")
#index_array=("min_fidelity_encode_ratio") 


for i in ${index_array[@]}
do
    #rm $i/*eps
    python3 $i/ChartGenerator.py
    cp $i/*.eps eps/QEC_$i.eps 
    cp $i/*.jpg jpg/QEC_$i.jpg 
    
done