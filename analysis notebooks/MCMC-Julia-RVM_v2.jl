#import Pkg 
#Pkg.add("AdaptiveMCMC")
#Pkg.add("Interpolations")
#Pkg.add("FLoops")
#Pkg.add(url="https://github.com/mvihola/AdaptiveMCMC.jl")
#Pkg.instantiate()

print("Running on ",Threads.nthreads()," threads\n")

using DelimitedFiles, Interpolations, AdaptiveMCMC


workdir="./"
juliadir="./"
foro = readdlm(workdir*"photondata.txt")

p0= [8.27378357e-02, 50, 75, 10, 0.5]
niter=1_000
nstep=2_000

include(juliadir*"/MaxLike_RVM.jl")

Sp=nothing
Rp=nothing
for i in 1:nstep
    out4 = adaptive_rwm(p0, 
                    x->MaxLike_RVM(x,foro), 
                    niter; b=0, algorithm=:am, Sp=Sp, Rp=Rp)

    allout=out4.X
    global Sp=out4.S
    global Rp=out4.R

    allout=out4.X[:,2:end]
    global p0=allout[:,end]
    print(i*niter,":",p0,"\n")
    print(Sp)
    print(Rp)
    open("minimum2M_RVM_v2.txt", "a") do io
           writedlm(io, transpose(allout))
    end
end
