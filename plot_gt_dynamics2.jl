function plot_gt_dynamics2(dishes,gt_data,trait_one,trait_two,p_ext_mix,axes_name)

    end_one = zeros(5)
    end_two = zeros(5)

    for i = 1:5
        dish_no = dishes[i]
        segment_one = trait_one[findall(gt_data.Dish_Number .== dish_no)]
        segment_two = trait_two[findall(gt_data.Dish_Number .== dish_no)]
        segment_days = days[findall(gt_data.Dish_Number .== dish_no)]
        seg_df = DataFrame(s1 = segment_one[:],s2 = segment_two[:],s3 = segment_days[:])
        sort!(seg_df,:s3)
        end_one[i] = segment_one[end]
        end_two[i] = segment_two[end] 
    end
    
    scatter!(axes_name,end_one[findall(p_ext_mix .< 18)],end_two[findall(p_ext_mix .< 18)],color=:black)
    scatter!(axes_name,end_one[findall(p_ext_mix .> 17)],end_two[findall(p_ext_mix .> 17)],color=:yellow)

return end_one, end_two
end