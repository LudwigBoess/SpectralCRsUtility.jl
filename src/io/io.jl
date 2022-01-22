using GadgetUnits
using GadgetIO
using Base.Threads

function readSingleCRShockDataFromOutputFile(file::String)

        # read file into memory
        f = open(file)
        lines = readlines(f)
        close(f)

        # filter only relevant lines. Output of every line starts with "CR DATA".
        g = occursin.("CR DATA", lines)
        lines = lines[g]

        # get number of output lines
        N = length(lines)

        # init datatype for analyze
        cr = CRShockData(N)

        for i ∈ 1:N
                line_content            = split(lines[i])
                cr.dt[i]                = parse(Float64,line_content[3])
                cr.Mach[i]              = parse(Float64,line_content[4])
                cr.Shock_Speed[i]       = parse(Float64,line_content[5])
                cr.Shock_Compress[i]    = parse(Float64,line_content[6])
                cr.Shock_Energy_In[i]   = parse(Float64,line_content[7])
                cr.Shock_Energy_Real[i] = parse(Float64,line_content[8])
                cr.Energy_P[i]          = parse(Float64,line_content[9])
                cr.Energy_e[i]          = parse(Float64,line_content[10])
        end

        return cr
end



"""
    getCRMomentumDistributionFromPartID( snap_file::String, ID::Integer;
                                         pmin::Real=1.0, pmax::Real=1.0e6,
                                         Nbins::Integer=0, mode::Int64=3)

Reads the spectra from a single SPH particle via the particle ID.
"""
function getCRMomentumDistributionFromPartID(snap_file::String, ID::Integer;
                                             pmin::Real=1.0, pmax::Real=1.0e6,
                                             Nbins::Integer=0, mode::Int64=3,
                                             protons::Bool=true, electrons::Bool=true)

    h = head_to_obj(snap_file)
    info = read_info(snap_file)

    if info == 1
        if Nbins == 0
            error("Can't read spectrum! No info block present!\nSupply number of momentum bins to proceed!")
        else
            info = Array{InfoLine,1}(undef,7)
            info[1] = InfoLine("ID",    UInt32, Int32(1),     [1, 0, 0, 0, 0, 0])
            info[2] = InfoLine("CRpN", Float32, Int32(Nbins), [1, 0, 0, 0, 0, 0])
            info[3] = InfoLine("CRpS", Float32, Int32(Nbins), [1, 0, 0, 0, 0, 0])
            info[4] = InfoLine("CRpC", Float32, Int32(1),     [1, 0, 0, 0, 0, 0])
            info[5] = InfoLine("CReN", Float32, Int32(Nbins), [1, 0, 0, 0, 0, 0])
            info[6] = InfoLine("CReS", Float32, Int32(Nbins), [1, 0, 0, 0, 0, 0])
            info[7] = InfoLine("CReC", Float32, Int32(1),     [1, 0, 0, 0, 0, 0])
        end
    end

    # read block positions to speed up IO
    block_positions = GadgetIO.get_block_positions(snap_file)

    id = read_block(snap_file, "ID",
                    info=info[getfield.(info, :block_name) .== "ID"][1],
                    parttype=0, block_position=block_positions["ID"])

    # select the position of the requested ID
    part = findfirst( id .== UInt32(ID) )[1]

    # protons
    if protons
        CRpN = Float64.(
            read_block(snap_file, "CRpN",
                        info=info[getfield.(info, :block_name) .== "CRpN"][1],
                        parttype=0, block_position=block_positions["CRpN"])[:,part]
                )
        CRpS = read_block(snap_file, "CRpS",
                        info=info[getfield.(info, :block_name) .== "CRpS"][1],
                        parttype=0, block_position=block_positions["CRpS"])[:,part]
        CRpC = read_block(snap_file, "CRpC",
                        info=info[getfield.(info, :block_name) .== "CRpC"][1],
                        parttype=0, block_position=block_positions["CRpC"])[part]

        Nbins = size(CRpS,1)
        par = CRMomentumDistributionConfig(pmin, pmax, Nbins, mode)
        CRp = CRMomentumDistribution( CRpN, CRpS, CRpC, par.pmin, par.pmax, par.mc_e )

    end

    # electrons
    if electrons
        CReN = Float64.(
            read_block(snap_file, "CReN",
                        info=info[getfield.(info, :block_name) .== "CReN"][1],
                        parttype=0, block_position=block_positions["CReN"])[:,part]
                )
        CReS = read_block(snap_file, "CReS",
                        info=info[getfield.(info, :block_name) .== "CReS"][1],
                        parttype=0, block_position=block_positions["CReS"])[:,part]
        CReC = read_block(snap_file, "CReC",
                        info=info[getfield.(info, :block_name) .== "CReC"][1],
                        parttype=0, block_position=block_positions["CReC"])[part]

        Nbins = size(CReS,1)
        par = CRMomentumDistributionConfig(pmin, pmax, Nbins, mode)
        CRe = CRMomentumDistribution( CReN, CReS, CReC, par.pmin, par.pmax, par.mc_p )
    end

    if protons && electrons
        return CRp, CRe      
    elseif protons 
        return CRp
    elseif electrons
        return CRe
    end
end



"""
    get_synchrotron_spectrum_from_part_id(snap_file::String, ID::Integer, ν_array;
                                          pmin::Real=1.0, pmax::Real=1.0e6,
                                          Nbins::Integer=0)

Reads the synchrotron spectrum from a single SPH particle via the particle ID.
"""
function get_synchrotron_spectrum_from_part_id(snap_file::String, ID::Integer, ν_array, B = 0.0;
                                             pmin::Real=1.0, pmax::Real=1.0e6,
                                             Nbins::Integer=0)

    h = head_to_obj(snap_file)
    info = read_info(snap_file)

    if info == 1
        if Nbins == 0
            error("Can't read spectrum! No info block present!\nSupply number of momentum bins to proceed!")
        else
            info = Array{InfoLine,1}(undef,7)
            info[1] = InfoLine("ID",    UInt32, Int32(1),     [1, 0, 0, 0, 0, 0])
            info[5] = InfoLine("CReN", Float32, Int32(Nbins), [1, 0, 0, 0, 0, 0])
            info[6] = InfoLine("CReS", Float32, Int32(Nbins), [1, 0, 0, 0, 0, 0])
            info[7] = InfoLine("CReC", Float32, Int32(1),     [1, 0, 0, 0, 0, 0])
        end
    end

    # read block positions to speed up IO
    block_positions = GadgetIO.get_block_positions(snap_file)

    id = read_block(snap_file, "ID",
                    info=info[getfield.(info, :block_name) .== "ID"][1],
                    parttype=0, block_position=block_positions["ID"])

    # select the position of the requested ID
    part = findfirst( id .== UInt32(ID) )[1]

    # define unit struct 
    GU = GadgetPhysical(h)

    # read CR data
    CReN = GU.CR_norm .* 10.0.^read_block(snap_file, "CReN",
                                          info=info[getfield.(info, :block_name) .== "CReN"][1],
                                          parttype=0, block_position=block_positions["CReN"])[:,part]
            
    CReS = read_block(snap_file, "CReS",
                    info=info[getfield.(info, :block_name) .== "CReS"][1],
                    parttype=0, block_position=block_positions["CReS"])[:,part] .|> Float64
    CReC = read_block(snap_file, "CReC",
                    info=info[getfield.(info, :block_name) .== "CReC"][1],
                    parttype=0, block_position=block_positions["CReC"])[part] .|> Float64

    Nbins = size(CReS,1)

    # read magnetic field (if present)
    if B == 0.0
        bfld  = read_block(snap_file, "BFLD", parttype=0)[:,part]
        B = √( bfld[1]^2 + bfld[2]^2 + bfld[3]^2 )
    end
    
    par = CRMomentumDistributionConfig(pmin, pmax, Nbins)

    j_ν = Vector{Float64}(undef, length(ν_array))
    @sync begin 
        @inbounds for i = 1:length(ν_array)
            @spawn begin 
                j_ν[i] = synchrotron_emission(CReN, CReS, CReC, B, par, ν0 = ν_array[i], 
                                        reduce_spectrum=true,
                                        integrate_pitch_angle = true)
            end
        end
    end

    return j_ν
    
end




"""
    write_crp_cre_to_txt( t::Vector{<:Real}, CRp::CRMomentumDistribution, CRe::CRMomentumDistribution, 
                          output_file::String )

Write CR Proton and Electron spectra for a series of time steps `t` to a txt file.
"""
function write_crp_cre_to_txt(t::Vector{<:Real}, CRp::CRMomentumDistribution, CRe::CRMomentumDistribution, 
                              output_file::String)

    data = Matrix{Float64}(undef, length(t), 1+4*length(CRp[1].norm))

    for i = 1:length(t)
        data[i,:] = [t[i] CRp[i].bound[1:end-1]' CRp[i].norm' CRe[i].bound[1:end-1]' CRe[i].norm' ]
    end

    writedlm(output_file, data)
end

"""
    read_crp_cre_from_txt(filename::String)

Read CR Proton and Electron spectra from a txt file.
"""
function read_crp_cre_from_txt(filename::String)

    data = readdlm(filename)

    N = size(data,1)
    Nbins = Int64((size(data,2) - 1) / 8)

    CRp = Array{CRMomentumDistribution,1}(undef,N)
    CRe = Array{CRMomentumDistribution,1}(undef,N)

    t = data[:,1]

    for i = 1:N
        CRp[i] = CRMomentumDistribution([data[i,2:2Nbins+1]; data[i,2Nbins+1]], data[i,2Nbins+2:4Nbins+1])
        CRe[i] = CRMomentumDistribution([data[i,4Nbins+2:6Nbins+1]; data[i,6Nbins+1]], data[i,6Nbins+2:8Nbins+1])
    end

    return t, CRp, CRe
end


function write_cr_to_txt(t::Vector{<:Real}, CR::Vector{CRMomentumDistribution}, output_file::String)

    data = Matrix{Float64}(undef, length(t), 1+2*length(CR[1].norm))

    for i = 1:length(t)
        data[i,:] = [t[i] CR[i].bound[1:end-1]' CR[i].norm' ]
    end

    writedlm(output_file, data)
end

function read_cr_from_txt(fi)

    data = readdlm(fi)

    N = size(data,1)
    Nbins = Int64((size(data,2) - 1) / 4)

    CR = Array{CRMomentumDistribution,1}(undef,N)

    t = data[:,1]

    for i = 1:N
        CR[i] = CRMomentumDistribution([data[i,2:2Nbins+1]; data[i,2Nbins+1]], data[i,2Nbins+2:4Nbins+1])
    end

    return t, CR
end