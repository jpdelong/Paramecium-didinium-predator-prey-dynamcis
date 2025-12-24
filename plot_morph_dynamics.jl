function plot_morph_dynamics(num_reps,morphdata,mix_dishes,p_ext_mix,trait_one,trait_two,axes_name,
    GT_starts,find_col1,find_col2)

    mix_end_one = zeros(1,num_reps)
    mix_end_two = zeros(1,num_reps)

    for i = 1:num_reps
        dish_no = mix_dishes[i]
        segment_one = trait_one[findall(morphdata.Dish_Number .== dish_no)]
        segment_two = trait_two[findall(morphdata.Dish_Number .== dish_no)]
        lines!(axes_name,segment_one,segment_two,color=:gray)
        scatter!(axes_name,segment_one[1],segment_two[1],color=col_mix,
           marker=:star5,markersize=20)
        mix_end_one[i] = segment_one[end]
        mix_end_two[i] = segment_two[end]
    end
    scatter!(axes_name,mix_end_one[findall(p_ext_mix .< 22)],mix_end_two[findall(p_ext_mix .< 22)],color=:black)
    scatter!(axes_name,mix_end_one[findall(p_ext_mix .> 21)],mix_end_two[findall(p_ext_mix .> 21)],color=:white,strokecolor=:black,strokewidth=1)

    #scatter!(axes_name,mix_end_one,mix_end_two,color=p_ext_mix,colormap=:roma)
    #scatter!(axes_name,mix_end_one,mix_end_two,color=:gray)
    #scatter!(axes_name,mix_end_one,mix_end_two,color=:gray)

    #GT_start_trait1 = GT_starts[:,find_col1]
    #GT_start_trait2 = GT_starts[:,find_col2]

    #scatter!(axes_name,GT_start_trait1,GT_start_trait2,
        #color=[col_9,col_33,col_34,col_56,col_57],marker=:rect,markersize=20)
    
    #scatter!(axes_name,g9_data[:,find_col1],g9_data[:,find_col2],color=col_9)
    #scatter!(axes_name,g33_data[:,find_col1],g33_data[:,find_col2],color=col_33)

return mix_end_one, mix_end_two

end
