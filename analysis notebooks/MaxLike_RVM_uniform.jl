using FLoops

#=
    param:
     1 Q/I            - param[1] (non-negative)
     2 alphap         - line of sight and rotation axis
     3 betad          - colatitude of dipole axis
     4 posang         - position angle of rotation axis
     5 phid           - longitude of dipole axis (0 to 1)

    data: matrix containing phase, evtlist['Q']*modf, evtlist['U']*modf
=#

function MaxLike_RVM(param::Array{Float64}, data::Matrix{Float64})
    if param[1]<0 || param[1]>1 || param[2]<0 || param[2]>180 || param[3]>90 || param[3]<0 || param[4]<0 || param[4]>180 || param[5]<0 || param[5]>1        
        return -Inf
    end
    l=size(data)[2]
    radang=param[4]/180*pi
    polfrac=param[1]
    halfamb=(param[2]-param[3])/2/180*pi
    halfapb=(param[2]+param[3])/2/180*pi
    shalfamb=sin(halfamb)
    chalfamb=cos(halfamb)
    shalfapb=sin(halfapb)
    chalfapb=cos(halfapb)
    @floop for i in 1:l
        @inbounds tanhalfC=tan((data[1,i]-param[5])*pi)
        halfAmB=atan(shalfamb,shalfapb*tanhalfC)
        halfApB=atan(chalfamb,chalfapb*tanhalfC)
        twoang=2*((halfApB-halfAmB)+radang)
        qloc=cos(twoang)
        uloc=sin(twoang)
        @inbounds @reduce ( c+=log(1+polfrac*(qloc*data[2,i]+uloc*data[3,i])) )
    end 
    return c
end
