function plot_gt_dynamics(dishes,gt_data,trait_one,trait_two,axes_name,gt_color,extinction_day)

    end_one = zeros(5)
    end_two = zeros(5)

    for i = 1:5
        dish_no = dishes[i]
        segment_one = trait_one[findall(gt_data.Dish_Number .== dish_no .&& gt_data.Day_of .< extinction_day)]
        segment_two = trait_two[findall(gt_data.Dish_Number .== dish_no .&& gt_data.Day_of .< extinction_day)]
        segment_days = gt_data.Day_of[findall(gt_data.Dish_Number .== dish_no .&& gt_data.Day_of .< extinction_day)]
        seg_df = DataFrame(s1 = segment_one[:],s2 = segment_two[:],s3 = segment_days[:])
        sort!(seg_df,:s3)

        lines!(axes_name,seg_df.s1,seg_df.s2,color=gt_color,alpha=0.5)
        scatter!(axes_name,segment_one[1],segment_two[1],color=gt_color,
            marker=:rect,markersize=20)
        end_one[i] = segment_one[end]
        end_two[i] = segment_two[end] 
    end

return end_one, end_two
end