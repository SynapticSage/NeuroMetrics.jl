jerryriggedsuperanimal = false
if dojerryriggedsuperanimal
    bins=61
    k=:all
    Plot.setfolder("jerry-rigged superanimal histogram ")
    H=[]
    @showprogress for k in (:all, :nontask, :cue, :mem_error, :cue_error, :task,:memory, :mem_correct, :cue_correct)
        key1, key2 = keyz1[k], keyz2[k]
        k  = string(k)
        f1,f2   = prep_ry16(matrixform(F[key1])), prep_ry22(matrixform(F[key2]))
        @info "$key" #"$(length(unique(f[:unit])))"

        vals1,vals2 = f1[shift=At(0)][:bestshift_bitsperspike], f2[shift=At(0)][:bestshift_bitsperspike]
        edges = LinRange(minimum(f[:shift]), maximum(f[:shift]), bins+1)
        xlim = (edges[begin], edges[end])

        #bl = bayesian_blocks(vals, prior=AIC(), resolution=1000)
        #support, density = to_pdf(bl)
        #plot(support, density)

       h= histogram(vcat(vals1, vals2); bins, xlim, title=k, label="")
       push!(H,h)
        vline!([0],c=:black,linestyle=:dash)
        Plot.save(k)
    end

    bins=61
    k=:all
    Plot.setfolder("jerry-rigged superanimal histogram  comparison")
    Plot.setappend("animal=jerryrig")
    H=[]
    @showprogress for (k1,k2) in ((:cue,), (:memory,))

        key11, key12 = keyz1[k1], keyz2[k1]
        key21, key22 = keyz1[k2], keyz2[k2]

        k1, k2  = string(k1),string(k2)
        f11 = prep_ry16(matrixform(F[key11]))
        f12 = prep_ry22(matrixform(F[key12]))
        f21 = prep_ry16(matrixform(F[key21]))
        f22 = prep_ry22(matrixform(F[key22]))

        #@info "$key" #"$(length(unique(f[:unit])))"

        vals11,vals12,vals21,vals22 = f11[shift=At(0)][:bestshift_bitsperspike], 
                        f12[shift=At(0)][:bestshift_bitsperspike],
                        f21[shift=At(0)][:bestshift_bitsperspike],
                        f22[shift=At(0)][:bestshift_bitsperspike]

        edges = LinRange(minimum(f[:shift]), maximum(f[:shift]), bins+1)
        xlim = (edges[begin], edges[end])

       #bl = bayesian_blocks(vals, prior=AIC(), resolution=1000)
       #support, density = to_pdf(bl)
       #plot(support, density)

      V1=vcat(vals11, vals12)
      V2=vcat(vals21, vals22)
      h1 = histogram(V1;  bins=9, xlim, title=k1 * " " * k2, label=k1, alpha=0.5)
      h2 = histogram!(V2; bins=9, xlim, label=k2, alpha=0.5)
      Plot.save((;k1,k2))
      using HypothesisTests
      e1, e2 = ecdf(V1), ecdf(V2)
      K = KruskalWallisTest(e1.sorted_values, e2.sorted_values)
      plot(e1,label=k1)
      plot!(e2, title="pval=$(pvalue(K))", xlabel="best shift on bits per spike", label=k2)
      Plot.save("bits per spike BEST SHIFT IS NOT sig different")


      push!(H,h)
       vline!([0],c=:black,linestyle=:dash)
       Plot.save(k)
    end
end
