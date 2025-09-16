function end_points(axes_name,end_one,end_two,p_ext,cutoff)

scatter!(axes_name,end_one[findall(p_ext .< cutoff+1)],end_two[findall(p_ext .< cutoff+1)],color=:black)
scatter!(axes_name,end_one[findall(p_ext .> cutoff)],end_two[findall(p_ext .> cutoff)],color=:yellow)

end