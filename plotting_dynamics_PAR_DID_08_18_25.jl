dir = "C:/Users/jdelong2/OneDrive - University of Nebraska/Projects in progress/Paramecium genotypes Didinium"
cd(dir)

##
using DataFrames, CSV, Statistics, MultivariateStats, GLM
using GLMakie
using Tables
using Biplots

include("stats_on_reps.jl")

# the structure of this code is to open up figures, create a layout, and
# then plot stuff iteratively into the layouts
# all the figure structure is done in advance
# then for each treatment, we grab some data and drop it on the appropriate 
# "frame"

# start with a choice of color schemes
gradient = cgrad(:darktest,6,categorical = true) # create a color gradient
col_mix = gradient[1]
col_9 = gradient[2]
col_33 = gradient[3]
col_34 = gradient[4]
col_56 = gradient[5]
col_57 = gradient[6]

# open up the frames and name the subpanels with a call to Axis
# Axis specifies the plot index location [row,col]
# we can add some layout details in here too like labels

f = Figure(size = (1000, 400))
    ax_56p = Axis(f[1, 1], title = "Paramecium")
    ax_56d = Axis(f[1, 2], title = "Didinium")
    ax_57p = Axis(f[2, 1])
    ax_57d = Axis(f[2, 2])
    ax_34p = Axis(f[3, 1], ylabel = "Cells per mL")
    ax_34d = Axis(f[3, 2])
    ax_9p = Axis(f[4, 1])
    ax_9d = Axis(f[4, 2])
    ax_33p = Axis(f[5, 1])
    ax_33d = Axis(f[5, 2])
    ax_para = Axis(f[6, 1], xlabel = "Days")
    ax_didi = Axis(f[6, 2], xlabel = "Days")

    # this lines up the ticks across subplanels
    linkxaxes!(ax_56p, ax_57p, ax_34p, ax_9p, ax_33p)
    linkxaxes!(ax_56d, ax_57d, ax_34d, ax_9d, ax_33d)

f2 = Figure()
    ax_SS = Axis(f2[1, 1], ylabel = "Paramecium density", ylabelsize=20,
        xlabel = "Didinium density", xlabelsize=20)
    
f3 = Figure()
    ax_days_ext = Axis(f3[1, 1], ylabel = "Day of extinction",ylabelsize=20,
        xticks=([1,2,3,4,5,6], ["Mixed","9","33","34","56","57"]))
    ax_cv = Axis(f3[2, 1], ylabel = "CV of density", xlabel = "Genotype treatment (line)",
        xlabelsize=20,ylabelsize=20,xticks=([1,2,3,4,5,6], ["Mixed","9","33","34","56","57"]))

# set up x axis positions by treatment for the bar and violin plots
p_pos = [0.9,1.9,2.9,3.9,4.9,5.9]
d_pos = [1.1,2.1,3.1,4.1,5.1,6.1]


#####
# start with GT 56
pddata = CSV.read("Did_par_56_all.csv",DataFrame)
times = pddata[:,:Day] .+ 0.0
Didi_density = pddata[:,2:6]
Para_density = pddata[:,7:11]
Para_average = pddata[:,12]
Didi_average = pddata[:,13]

    series!(ax_56p,times, Matrix(Para_density)', solid_color=col_56, linewidth=1)
    c1 = series!(ax_56p,times,Para_average', solid_color=col_56,label="Line 56",
        linewidth=4)
    axislegend(ax_56p,[c1],["Line 56"], position = :rt)

    series!(ax_56d,times, Matrix(Didi_density)', solid_color=col_56, linewidth=1)
    series!(ax_56d,times,Didi_average', solid_color=col_56,label="",
        linewidth=4)
    z1 = series!(ax_SS,Didi_average,Para_average', solid_color=col_56,label="Line 56",
        linewidth=3)

    scatter!(ax_SS,Didi_average[1],Para_average[1],color=:black,
        markersize = 20)

    # call this function for some rep information
    rep_results = stats_on_reps(Didi_density,Para_density,times,5)
    day_ext_56 = rep_results[:,3]
    d_num_ext = sum(rep_results[:,1] .< rep_results[:,2]) # count how many reps in which D went extinct first
    p_num_ext = sum(rep_results[:,1] .> rep_results[:,2]) # count how many reps in which P went extinct first
    p_dens_cv_56 = rep_results[:,4]
    d_dens_cv_56 = rep_results[:,5]

    # now plot the days to extinction
    plot!(ax_days_ext,p_pos[5].*ones(p_num_ext) .+ 0.1 .* rand.(),day_ext_56[rep_results[:,1] .> rep_results[:,2]],color=col_56,marker=:ltriangle)
    plot!(ax_days_ext,d_pos[5].*ones(d_num_ext) .+ 0.1 .* rand.(),day_ext_56[rep_results[:,1] .< rep_results[:,2]],color=col_56)

    # now cv of abundances
    plot!(ax_cv,p_pos[5].*ones(5) .+ 0.1 .* rand.(),p_dens_cv_56,color=col_56,marker=:ltriangle)
    plot!(ax_cv,d_pos[5].*ones(5) .+ 0.1 .* rand.(),d_dens_cv_56,color=col_56)

#####
# GT 57
pddata = CSV.read("Did_par_57_all.csv",DataFrame)
times = pddata[:,:Day]
Didi_density = pddata[:,2:6]
Para_density = pddata[:,7:11]
Para_average = pddata[:,12]
Didi_average = pddata[:,13]

    series!(ax_57p,times, Matrix(Para_density)', solid_color=col_57, linewidth=1)
    c2 = series!(ax_57p,times,Para_average', solid_color=col_57,label="Line 57",
        linewidth=4)
    axislegend(ax_57p,[c2],["Line 57"], position = :rt)

    series!(ax_57d,times, Matrix(Didi_density)', solid_color=col_57, linewidth=1)
    series!(ax_57d,times,Didi_average', solid_color=col_57,label="",
        linewidth=4)

    z2 = series!(ax_SS,Didi_average,Para_average', solid_color=col_57,label="Line 57",
        linewidth=3)

    # call this function for some rep information
    rep_results = stats_on_reps(Didi_density,Para_density,times,5)
    day_ext_57 = rep_results[:,3]
    d_num_ext = sum(rep_results[:,1] .< rep_results[:,2]) # count how many reps in which D went extinct first
    p_num_ext = sum(rep_results[:,1] .> rep_results[:,2]) # count how many reps in which P went extinct first
    p_dens_cv_57 = rep_results[:,4]
    d_dens_cv_57 = rep_results[:,5]
    
    # now plot the days to extinction
    plot!(ax_days_ext,p_pos[6].*ones(p_num_ext) .+ 0.1 .* rand.(),day_ext_57[rep_results[:,1] .> rep_results[:,2]],color=col_57,marker=:ltriangle)
    plot!(ax_days_ext,d_pos[6].*ones(d_num_ext) .+ 0.1 .* rand.(),day_ext_57[rep_results[:,1] .< rep_results[:,2]],color=col_57)

    # now cv of abundances
    plot!(ax_cv,p_pos[6].*ones(5) .+ 0.1 .* rand.(),p_dens_cv_57,color=col_57,marker=:ltriangle)
    plot!(ax_cv,d_pos[6].*ones(5) .+ 0.1 .* rand.(),d_dens_cv_57,color=col_57)

#####
# GT 34
pddata = CSV.read("Did_par_34_all.csv",DataFrame)
times = pddata[:,:Day]
Didi_density = pddata[:,2:6]
Para_density = pddata[:,7:11]
Para_average = pddata[:,12]
Didi_average = pddata[:,13]

    series!(ax_34p,times, Matrix(Para_density)', solid_color=col_34, linewidth=1)
    c3 = series!(ax_34p,times,Para_average', solid_color=col_34,label="Line 34",
        linewidth=4)
    axislegend(ax_34p,[c3],["Line 34"], position = :rt)

    series!(ax_34d,times, Matrix(Didi_density)', solid_color=col_34, linewidth=1)
    series!(ax_34d,times,Didi_average', solid_color=col_34,label="",
        linewidth=4)

    z3 = series!(ax_SS,Didi_average,Para_average', solid_color=col_34,label="Line 34",
        linewidth=3)

    # call this function for some rep information
    rep_results = stats_on_reps(Didi_density,Para_density,times,5)
    day_ext_34 = rep_results[:,1]
    d_num_ext = sum(rep_results[:,1] .< rep_results[:,2]) # count how many reps in which D went extinct first
    p_num_ext = sum(rep_results[:,1] .> rep_results[:,2]) # count how many reps in which P went extinct first
    p_dens_cv_34 = rep_results[:,4]
    d_dens_cv_34 = rep_results[:,5]

    # now plot the days to extinction
    plot!(ax_days_ext,p_pos[4].*ones(p_num_ext) .+ 0.1 .* rand.(),day_ext_34[rep_results[:,1] .> rep_results[:,2]],color=col_34,marker=:ltriangle)
    plot!(ax_days_ext,d_pos[4].*ones(d_num_ext) .+ 0.1 .* rand.(),day_ext_34[rep_results[:,1] .< rep_results[:,2]],color=col_34)

    # now cv of abundances
    plot!(ax_cv,p_pos[4].*ones(5) .+ 0.1 .* rand.(),p_dens_cv_34,color=col_34,marker=:ltriangle)
    plot!(ax_cv,d_pos[4].*ones(5) .+ 0.1 .* rand.(),d_dens_cv_34,color=col_34)


# GT 9
pddata = CSV.read("Did_par_9_all.csv",DataFrame)
times = pddata[:,:Day]
Didi_density = pddata[:,2:6]
Para_density = pddata[:,7:11]
Para_average = pddata[:,12]
Didi_average = pddata[:,13]

    series!(ax_9p,times, Matrix(Para_density)', solid_color=col_9, linewidth=1)
    c4 = series!(ax_9p,times,Para_average', solid_color=col_9,label="Line 9",
        linewidth=4)
    axislegend(ax_9p,[c4],["Line 9"], position = :rt)

    series!(ax_9d,times, Matrix(Didi_density)', solid_color=col_9, linewidth=1)
    series!(ax_9d,times,Didi_average', solid_color=col_9,label="",
        linewidth=4)

    z4 = series!(ax_SS,Didi_average,Para_average', solid_color=col_9,label="",
        linewidth=3)

    # call this function for some rep information
    rep_results = stats_on_reps(Didi_density,Para_density,times,5)
    day_ext_9 = rep_results[:,1]
    d_num_ext = sum(rep_results[:,1] .< rep_results[:,2]) # count how many reps in which D went extinct first
    p_num_ext = sum(rep_results[:,1] .> rep_results[:,2]) # count how many reps in which P went extinct first
    p_dens_cv_9 = rep_results[:,4]
    d_dens_cv_9 = rep_results[:,5]

    # now plot the days to extinction
    plot!(ax_days_ext,p_pos[2].*ones(p_num_ext) .+ 0.1 .* rand.(),day_ext_9[rep_results[:,1] .> rep_results[:,2]],color=col_9,marker=:ltriangle)
    plot!(ax_days_ext,d_pos[2].*ones(d_num_ext) .+ 0.1 .* rand.(),day_ext_9[rep_results[:,1] .< rep_results[:,2]],color=col_9)

    # now cv of abundances
    plot!(ax_cv,p_pos[2].*ones(5) .+ 0.1 .* rand.(),p_dens_cv_9,color=col_9,marker=:ltriangle)
    plot!(ax_cv,d_pos[2].*ones(5) .+ 0.1 .* rand.(),d_dens_cv_9,color=col_9)


#####
# GT 33
pddata = CSV.read("Did_par_33_all.csv",DataFrame)
times = pddata[:,:Day]
Didi_density = pddata[:,2:6]
Para_density = pddata[:,7:11]
Para_average = pddata[:,12]
Didi_average = pddata[:,13]

    series!(ax_33p,times, Matrix(Para_density)', solid_color=col_33, linewidth=1)
    c5 = series!(ax_33p,times,Para_average', solid_color=col_33,label="Line 33",
        linewidth=4)
    axislegend(ax_33p,[c5],["Line 33"], position = :rt)

    series!(ax_33d,times, Matrix(Didi_density)', solid_color=col_33, linewidth=1)
    series!(ax_33d,times,Didi_average', solid_color=col_33,label="",
        linewidth=4)

    z5 = series!(ax_SS,Didi_average,Para_average', solid_color=col_33,label="",
        linewidth=3)

    # call this function for some rep information
    rep_results = stats_on_reps(Didi_density,Para_density,times,5)
    day_ext_33 = rep_results[:,3] # find the time to extinction for either D or P
    d_num_ext = sum(rep_results[:,1] .< rep_results[:,2]) # count how many reps in which D went extinct first
    p_num_ext = sum(rep_results[:,1] .> rep_results[:,2]) # count how many reps in which P went extinct first
    p_dens_cv_33 = rep_results[:,4] # get cv for P
    d_dens_cv_33 = rep_results[:,5] # get cv for D

    # now plot the days to extinction
    plot!(ax_days_ext,p_pos[3].*ones(p_num_ext) .+ 0.1 .* rand.(),day_ext_33[rep_results[:,1] .> rep_results[:,2]],color=col_33,marker=:ltriangle)
    plot!(ax_days_ext,d_pos[3].*ones(d_num_ext) .+ 0.1 .* rand.(),day_ext_33[rep_results[:,1] .< rep_results[:,2]],color=col_33)

    # now cv of abundances
    plot!(ax_cv,p_pos[3].*ones(5) .+ 0.1 .* rand.(),p_dens_cv_33,color=col_33,marker=:ltriangle)
    plot!(ax_cv,d_pos[3].*ones(5) .+ 0.1 .* rand.(),d_dens_cv_33,color=col_33)

#####
# Mixed GT dishes
pddata = CSV.read("Mixed_genotype_dishes.csv",DataFrame)
times = pddata.Day
Didi_density = pddata[:,2:26]
Para_density = pddata[:,27:51]
Para_average = pddata[:,52]
Didi_average = pddata[:,53]

    series!(ax_para,times, Matrix(Para_density)', solid_color=col_mix, linewidth=1)
    c5 = series!(ax_para,times,Para_average', solid_color=col_mix,label="Mixed",
        linewidth=4)
    axislegend(ax_para,[c5],["Mixed"], position = :rt)

    series!(ax_didi,times, Matrix(Didi_density)', solid_color=col_mix, linewidth=1)
    series!(ax_didi,times,Didi_average', solid_color=col_mix,label="",
        linewidth=4)

    z6 = series!(ax_SS,Didi_average,Para_average', solid_color=col_mix,label="",
        linewidth=3)

    # call this function for some rep information
    rep_results = stats_on_reps(Didi_density,Para_density,times,25)
    day_ext_mix = rep_results[:,1]
    p_dens_cv_mix = rep_results[:,4]
    d_dens_cv_mix = rep_results[:,5]
    d_num_ext = sum(rep_results[:,1] .< rep_results[:,2]) # count how many reps in which D went extinct first
    p_num_ext = sum(rep_results[:,1] .> rep_results[:,2]) # count how many reps in which P went extinct first

    # now plot the days to extinction
    test1 = plot!(ax_days_ext,p_pos[1].*ones(p_num_ext) .+ 0.1 .* rand.(),day_ext_mix[rep_results[:,1] .> rep_results[:,2]],color=col_mix,marker=:ltriangle)
    test2 = plot!(ax_days_ext,d_pos[1].*ones(d_num_ext) .+ 0.1 .* rand.(),day_ext_mix[rep_results[:,1] .< rep_results[:,2]],color=col_mix)
    axislegend(ax_days_ext,[test1,test2],["P","D"], position = (0.25,0.5))

    # now cv of abundances
    test1 = plot!(ax_cv,p_pos[1].*ones(25) .+ 0.1 .* rand.(),p_dens_cv_mix,color=col_mix,marker=:ltriangle)
    test2 = plot!(ax_cv,d_pos[1].*ones(25) .+ 0.1 .* rand.(),d_dens_cv_mix,color=col_mix)

    axislegend(ax_SS,
    [z1, z2, z3, z4, z5, z6],
    ["GT_56", "GT_57", "GT_34", "GT_9", "GT_33", "Mixed"])

# get panel labels on plot f3
text!(ax_days_ext, (0.5, 21), text = "A", fontsize = 20)
text!(ax_cv, (0.5, 2), text = "B", fontsize = 20)


f
f2
f3

save("GT_timeseries.png",f)
save("GT_statespace.png",f2)
save("GT_stability.png",f3)

##############################################
# now do the trait dynamics plots
morphdata = CSV.read("MorphDataCombined.csv",DataFrame)

# include the plotting function
include("plot_morph_dynamics.jl")

# make a pca of all the traits in the morph data
morphdata_subset = Matrix(morphdata[:,9:32])
#names_subset = header_names([9:32])
#= for some reason have to write out the names in symbol format
names_subset=[:mean_area,:sd_area,:mean_perimeter,:sd_perimeter,:mean_major,
:sd_major,:mean_minor,:sd_minor,:mean_ar,:sd_ar,:sd_turning,:duration,:N_frames,
:max_net,:net_disp,:net_speed,:gross_disp,:gross_speed,:max_step,:min_step,:sd_step,:sd_gross_speed,:max_gross_speed,:min_gross_speed]=#
#:mean_turning ---- not using this because it has negative values - should be fixable?

names_subset=[:mean_area,:mean_major,:mean_minor,:mean_ar,:sd_turning,
    :net_disp,:gross_speed]

morphdata_subset = Matrix(morphdata[:, [:mean_area,:mean_major,:mean_minor,:mean_ar,:sd_turning,
    :net_disp,:gross_speed]])

    # from FUN ECOL paper, we just used: length, width, aspect ratio, mean turning angle, standard deviation of the turning angle, gross speed, net displacement, and standard deviation


M = fit(PCA, morphdata_subset'; maxoutdim=3)
pcadata = predict(M, morphdata_subset')'
pcadata = Matrix(pcadata)
# add the PCA axes to the dataframe
morphdata[!,:PCA1] .= pcadata[:,1]
morphdata[!,:PCA2] .= pcadata[:,2]

# or if we want the zscored PCA data
morphdata[!,:PCA1] .= (pcadata[:,1] .- mean(pcadata[:,1])) ./ std(pcadata[:,1])
morphdata[!,:PCA2] .= (pcadata[:,2] .- mean(pcadata[:,2])) ./ std(pcadata[:,2])

# 2D relative variation biplot with colored dots
table2 = (; zip(names_subset, eachcol(morphdata_subset))...)
fig2, ax2 = biplot(table2, kind = :rform, dotsize = 10, showdots=false)
# dotcolor = table2.mean_area, 
ax2.aspect = DataAspect()
save("GT_biplot.png",fig2)

# grab genotype-specific data and find start
g9_data = morphdata[occursin.("G9", morphdata.Genotype), :]
g9_start = g9_data[findfirst(==(1), g9_data.Day_of .== 0),:]
g33_data = morphdata[occursin.("G33", morphdata.Genotype), :]
g33_start = g33_data[findfirst(==(1), g33_data.Day_of .== 0),:]
g34_data = morphdata[occursin.("G34", morphdata.Genotype), :]
g34_start = g34_data[findfirst(==(1), g34_data.Day_of .== 0),:]
g56_data = morphdata[occursin.("G56", morphdata.Genotype), :]
g56_start = g56_data[findfirst(==(1), g56_data.Day_of .== 0),:]
g57_data = morphdata[occursin.("G57", morphdata.Genotype), :]
g57_start = g57_data[findfirst(==(1), g57_data.Day_of .== 0),:]

# pull the start rows together
GT_starts = append!(DataFrame(g9_start),DataFrame(g33_start),DataFrame(g34_start),
    DataFrame(g56_start),DataFrame(g57_start))
df_names = names(GT_starts)

# pull the dish id for the mixed genotype populations
find_mixes = morphdata[occursin.("Mix", morphdata.Genotype), :]
mix_dishes = unique(find_mixes.Dish_Number)
mix_dishes = mix_dishes[1:end]
num_reps = length(mix_dishes)

# open up a new figure
f4 = Figure(size = (1500, 500))
    ax_pca_10 = Axis(f4[1, 1],ylabel="PCA2",title = "Up to day 10")
    ax_pca_20 = Axis(f4[2, 1],ylabel="PCA2",title = "Up to day 20")
    ax_pca_25 = Axis(f4[3, 1],ylabel="PCA2",xlabel="PCA1",title = "Up to day 24")
    
    ax_pcade_10 = Axis(f4[1, 2],ylabel="Day of extinction",title = "Up to day 10")
    ax_pcade_20 = Axis(f4[2, 2],ylabel="Day of extinction",title = "Up to day 20")
    ax_pcade_25 = Axis(f4[3, 2],ylabel="Day of extinction",xlabel="PCA1",title = "Up to day 24")

# restrict data we use to FRs with more than a minimum amount of data
morphdata_10 = morphdata[findall(morphdata.Day_of .< 10),:]

    # create the plot for the PCA axes
    find_col1 = findfirst(==(1), df_names .== "PCA1")
    find_col2 = findfirst(==(1), df_names .== "PCA2")
    mix_end_one, mix_end_two = plot_morph_dynamics(num_reps,morphdata_10,mix_dishes,day_ext_mix,morphdata_10.PCA1,morphdata_10.PCA2,ax_pca_10,
        GT_starts,find_col1,find_col2)

    # make a dataset to test for an effect of PCA axes on extinction
    df_extinction = DataFrame(Days_to_extinction = day_ext_mix[:])
    df_extinction[!,:PCA1] .= mix_end_one[:]
    df_extinction[!,:PCA2] .= mix_end_two[:]

    lm_1 = lm(@formula(Days_to_extinction ~ PCA1 * PCA2), df_extinction)
    lm_1 = lm(@formula(Days_to_extinction ~ PCA1 + PCA2), df_extinction)
    lm_1 = lm(@formula(Days_to_extinction ~ PCA1), df_extinction)

    scatter!(ax_pcade_10,mix_end_one,day_ext_mix,color=:black)
    new_x = DataFrame(PCA1 = range(minimum(df_extinction.PCA1), maximum(df_extinction.PCA1), length=10))
    df_predict_y = predict(lm_1,new_x, interval=:confidence, level = 0.95)
    lines!(ax_pcade_10,new_x.PCA1, df_predict_y.prediction, linewidth = 2)
    band!(ax_pcade_10,new_x.PCA1,df_predict_y.lower,df_predict_y.upper; color = (:blue, 0.2))

# restrict data we use to FRs with more than a minimum amount of data
morphdata_20 = morphdata[findall(morphdata.Day_of .< 20),:]

    # create the plot for the PCA axes
    mix_end_one, mix_end_two = plot_morph_dynamics(num_reps,morphdata_20,mix_dishes,day_ext_mix,morphdata_20.PCA1,morphdata_20.PCA2,ax_pca_20,
        GT_starts,find_col1,find_col2)

    # make a dataset to test for an effect of PCA axes on extinction
    df_extinction = DataFrame(Days_to_extinction = day_ext_mix[:])
    df_extinction[!,:cv] .= p_dens_cv_mix[:]
    df_extinction[!,:PCA1] .= mix_end_one[:]
    df_extinction[!,:PCA2] .= mix_end_two[:]

    lm_1 = lm(@formula(Days_to_extinction ~ PCA1 * PCA2), df_extinction)
    lm_1 = lm(@formula(Days_to_extinction ~ PCA1 + PCA2), df_extinction)
    lm_1 = lm(@formula(Days_to_extinction ~ PCA1), df_extinction)

    scatter!(ax_pcade_20,mix_end_one,day_ext_mix,color=:black)
    new_x = DataFrame(PCA1 = range(minimum(df_extinction.PCA1), maximum(df_extinction.PCA1), length=10))
    df_predict_y = predict(lm_1,new_x, interval=:confidence, level = 0.95)
    lines!(ax_pcade_20,new_x.PCA1, df_predict_y.prediction, linewidth = 2)
    band!(ax_pcade_20,new_x.PCA1,df_predict_y.lower,df_predict_y.upper; color = (:blue, 0.2))

# no need to restrict data for the full time series

    # create the plot for the PCA axes
    mix_end_one, mix_end_two = plot_morph_dynamics(num_reps,morphdata,mix_dishes,day_ext_mix,morphdata.PCA1,morphdata.PCA2,ax_pca_25,
        GT_starts,find_col1,find_col2)

    # make a dataset to test for an effect of PCA axes on extinction
    df_extinction = DataFrame(Days_to_extinction = day_ext_mix[:])
    df_extinction[!,:cv] .= p_dens_cv_mix[:]
    df_extinction[!,:PCA1] .= mix_end_one[:]
    df_extinction[!,:PCA2] .= mix_end_two[:]

    lm_1 = lm(@formula(Days_to_extinction ~ PCA1 * PCA2), df_extinction)
    lm_1 = lm(@formula(Days_to_extinction ~ PCA1 + PCA2), df_extinction)
    lm_1 = lm(@formula(Days_to_extinction ~ PCA1), df_extinction)

    scatter!(ax_pcade_25,mix_end_one,day_ext_mix,color=:black)
    new_x = DataFrame(PCA1 = range(minimum(df_extinction.PCA1), maximum(df_extinction.PCA1), length=10))
    df_predict_y = predict(lm_1,new_x, interval=:confidence, level = 0.95)
    lines!(ax_pcade_25,new_x.PCA1, df_predict_y.prediction, linewidth = 2)
    band!(ax_pcade_25,new_x.PCA1,df_predict_y.lower,df_predict_y.upper; color = (:blue, 0.2))


    #=Colorbar(f4[1, 3], limits = (0, num_reps),
        colormap = cgrad(:roma, num_reps, categorical = true), size = num_reps)=#

    Legend(f4[1, 3], [PolyElement(color = col_9),PolyElement(color = col_33),PolyElement(color = col_34),
    PolyElement(color = col_56),PolyElement(color = col_57),PolyElement(color = :black),PolyElement(color = :yellow)],
    ["line 9","line 33","line 34","line 56","line 57","Extinct","Extant"])

    #["GT9-most extant","GT33-all extinct","GT34-half extinct","GT56-most extinct","GT57-most extinct","Extinct","Extant"])

f4
save("Stability_test.png",f4)



#= pull out variables of interest
days = morphdata[:,:Day_of] .+ 0.0
major_length = morphdata[:,:mean_major]
minor_length = morphdata[:,:mean_minor]
aspect_ratio = morphdata[:,:mean_ar]
mean_area = morphdata[:,:mean_area]
gross_speed = morphdata[:,:gross_speed]
mean_turning = morphdata[:,:mean_turning]
sd_turning = morphdata[:,:sd_turning]
max_gross_speed = morphdata[:,:max_gross_speed]
net_disp = morphdata[:,:net_disp]
dish = morphdata[:,:Dish_Number]
=#

##############################################
# now plot the phenotype dynamics for the individual GTs

# include the plotting function
include("plot_gt_dynamics.jl")
#include("plot_gt_dynamics2.jl")
include("end_points.jl")


f5 = Figure(size = (500, 500))
    g_length_ar = Axis(f5[1, 1],xlabel="PCA1",ylabel="PCA2")
    #end_length_ar = Axis(f5[2, 1],xlabel="PCA1",ylabel="PCA2")
        
    find_col1 = findfirst(==(1), df_names .== "PCA1")
    find_col2 = findfirst(==(1), df_names .== "PCA2")

    #=
    for i = 1:25
        dish_no = mix_dishes[i]
        segment_one = morphdata[findall(morphdata.Dish_Number .== dish_no),find_col1]
        segment_two = morphdata[findall(morphdata.Dish_Number .== dish_no),find_col2]
        lines!(g_length_ar,segment_one,segment_two,color=:black,alpha=1,linewidth=2)
    end=#

    extinction_day = 24.0

    ##### gt 9
    g9_dishes = unique(g9_data.Dish_Number)
    end_one, end_two = plot_gt_dynamics(g9_dishes,g9_data,g9_data[:,find_col1],g9_data[:,find_col2],g_length_ar,col_9,extinction_day)
    end_points(end_length_ar,end_one,end_two,day_ext_9,extinction_day)

    # make a dataset to test for an effect of PCA axes on extinction
    # got to build this after pulling data out for each GT in succession
    df_plasticity = DataFrame(PCA1 = end_one[:])
    df_plasticity[!,:PCA2] .= end_two[:]
    df_plasticity[!,:Day_extinct] .= day_ext_9[:]

    ##### gt 33
    g33_dishes = unique(g33_data.Dish_Number)
    end_one, end_two = plot_gt_dynamics(g33_dishes,g33_data,g33_data[:,find_col1],g33_data[:,find_col2],g_length_ar,col_33,extinction_day)
    end_points(end_length_ar,end_one,end_two,day_ext_33,extinction_day)

    # make a dataset to test for an effect of PCA axes on extinction
    df_plasticity2 = DataFrame(PCA1 = end_one[:])
    df_plasticity2[!,:PCA2] .= end_two[:]
    df_plasticity2[!,:Day_extinct] .= day_ext_33[:]

    g34_dishes = unique(g34_data.Dish_Number)
    end_one, end_two = plot_gt_dynamics(g34_dishes,g34_data,g34_data[:,find_col1],g34_data[:,find_col2],g_length_ar,col_34,extinction_day)
    end_points(end_length_ar,end_one,end_two,day_ext_34,extinction_day)

    # make a dataset to test for an effect of PCA axes on extinction
    df_plasticity3 = DataFrame(PCA1 = end_one[:])
    df_plasticity3[!,:PCA2] .= end_two[:]
    df_plasticity3[!,:Day_extinct] .= day_ext_34[:]

    g56_dishes = unique(g56_data.Dish_Number)
    end_one, end_two = plot_gt_dynamics(g56_dishes,g56_data,g56_data[:,find_col1],g56_data[:,find_col2],g_length_ar,col_56,extinction_day)
    end_points(end_length_ar,end_one,end_two,day_ext_56,extinction_day)

    # make a dataset to test for an effect of PCA axes on extinction
    df_plasticity4 = DataFrame(PCA1 = end_one[:])
    df_plasticity4[!,:PCA2] .= end_two[:]
    df_plasticity4[!,:Day_extinct] .= day_ext_56[:]

    g57_dishes = unique(g57_data.Dish_Number)
    end_one, end_two = plot_gt_dynamics(g57_dishes,g57_data,g57_data[:,find_col1],g57_data[:,find_col2],g_length_ar,col_57,extinction_day)
    end_points(end_length_ar,end_one,end_two,day_ext_57,extinction_day)

    # make a dataset to test for an effect of PCA axes on extinction
    df_plasticity5 = DataFrame(PCA1 = end_one[:])
    df_plasticity5[!,:PCA2] .= end_two[:]
    df_plasticity5[!,:Day_extinct] .= day_ext_57[:]

    Legend(f5[1, 2], [PolyElement(color = col_9),PolyElement(color = col_33),PolyElement(color = col_34),
    PolyElement(color = col_56),PolyElement(color = col_57)],
    ["line 9","line 33","line 34","line 56","line 57"])

f5
save("GT_phenotype_traj.png",f5)


    Legend(f5[2, 2], [PolyElement(color = :black),PolyElement(color = :yellow)],["Extinct","Extant"])


    df_plasticity_all = append!(df_plasticity,df_plasticity2,df_plasticity3,df_plasticity4,df_plasticity5)
    lm_1 = lm(@formula(Day_extinct ~ PCA1 * PCA2), df_plasticity_all)
    lm_1 = lm(@formula(Day_extinct ~ PCA1 + PCA2), df_plasticity_all)
