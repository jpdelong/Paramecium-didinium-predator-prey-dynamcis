function stats_on_reps(Didi_density,Para_density,times,reps)

p_ext = zeros(reps)
d_ext = zeros(reps)
p_dens_cv = zeros(reps)
d_dens_cv = zeros(reps)

for i = 1:reps
    # tally the day of extinction
    if isnothing(findfirst(==(0),skipmissing(Didi_density[:,i])))
        d_ext[i] = 24
    else    
        d_ext[i] = times[findfirst(==(0),skipmissing(Didi_density[:,i]))]
    end
    if isnothing(findfirst(==(0),skipmissing(Para_density[:,i])))
        p_ext[i] = 24
    else    
        p_ext[i] = times[findfirst(==(0),skipmissing(Para_density[:,i]))]
    end

    # calculate the cvs of both p and d abundances
    p_dens_cv[i] = std(skipmissing(Para_density[:,i])) / mean(skipmissing(Para_density[:,i]))
    d_dens_cv[i] = std(skipmissing(Didi_density[:,i])) / mean(skipmissing(Didi_density[:,i]))

end

return [d_ext p_ext min(d_ext,p_ext) p_dens_cv d_dens_cv]

end