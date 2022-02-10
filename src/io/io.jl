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



function get_particle_position_in_block(snap_file, ID)

    h = read_header(snap_file)
    
    Nfile = 0
    
    for i = 0:h.num_files-1

        filename = GadgetIO.select_file(snap_file, i)

        # read block positions to speed up IO
        id = read_block(filename, "ID", parttype=0)

        # select the position of the requested ID
        sel = findfirst( id .== UInt32(ID) )

        if !isnothing(sel)
            return i, sel[1]
        end
    end

    error("ID not found!")
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


    Nfile, sel = get_particle_position_in_block(snap_file, ID)

    # protons
    if protons

        filename = GadgetIO.select_file(snap_file, Nfile)

        CRpN = read_block(filename, "CRpN", parttype=0)[:,sel]
        CRpS = read_block(filename, "CRpS", parttype=0)[:,sel]
        CRpC = read_block(filename, "CRpC", parttype=0)[sel]

        Nbins = size(CRpS,1)
        par = CRMomentumDistributionConfig(pmin, pmax, Nbins, mode)
        CRp = CRMomentumDistribution( CRpN, CRpS, CRpC, par.pmin, par.pmax, par.mc_e )

    end

    # electrons
    if electrons
        
        filename = GadgetIO.select_file(snap_file, Nfile)

        CReN = read_block(filename, "CReN", parttype=0)[:,sel]
        CReS = read_block(filename, "CReS", parttype=0)[:,sel]
        CReC = read_block(filename, "CReC", parttype=0)[sel]

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
    getCRMomentumDistributionFromPartID( snap_file::String, ID::Integer;
                                         pmin::Real=1.0, pmax::Real=1.0e6,
                                         Nbins::Integer=0, mode::Int64=3)

Reads the spectra from a single SPH particle via the particle ID.
"""
function snapshot_to_spectra(snap_file::String;
                             pmin::Real=1.0, pmax::Real=1.0e6,
                             mode::Int64=3,
                             protons::Bool=true, electrons::Bool=true)

    h = read_header(snap_file)
    Npart =  h.nall[1]
    
    info = read_info(snap_file)

    if info == 1
        error("IO only works with INFO block.")
    end

    # protons
    if protons

        CRpP = read_block(snap_file, "CRpP", parttype=0)
        sel = findall(CRpP .> 0.0)
        Npart = length(sel)

        @info "Protons: $Npart / $(h.nall[1]) active"

        CRp  = Vector{CRMomentumDistribution}(undef, Npart)

        CRpN = read_block(snap_file, "CRpN", parttype=0)[:,sel]
        CRpS = read_block(snap_file, "CRpS", parttype=0)[:,sel]
        CRpC = read_block(snap_file, "CRpC", parttype=0)[sel]

        Nbins = size(CRpS,1)
        par = CRMomentumDistributionConfig(pmin, pmax, Nbins, mode)

        @sync begin
            @inbounds for i = 1:Npart
                @spawn begin
                    CRp[i] = CRMomentumDistribution( CRpN[:,i], CRpS[:,i], CRpC[i], par.pmin, par.pmax, par.mc_e )
                end
            end 
        end

    end

    # electrons
    if electrons

        CReP = read_block(snap_file, "CReP", parttype=0)
        sel = findall(CReP .> 0.0)
        Npart = length(sel)

        @info "Electrons: $Npart / $(h.nall[1]) active"

        CRe  = Vector{CRMomentumDistribution}(undef, Npart)

        CReN = read_block(snap_file, "CReN", parttype=0)[:,sel]
        CReS = read_block(snap_file, "CReS", parttype=0)[:,sel]
        CReC = read_block(snap_file, "CReC", parttype=0)[sel]

        Nbins = size(CRpS,1)
        par = CRMomentumDistributionConfig(pmin, pmax, Nbins, mode)

        @sync begin
            @inbounds for i = 1:Npart
                @spawn begin
                    CRe[i] = CRMomentumDistribution( CReN[:,i], CReS[:,i], CReC[i], par.pmin, par.pmax, par.mc_e )
                end
            end
        end
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
function get_synchrotron_spectrum_from_part_id(snap_file::String, ID::Integer, ν_array::Vector{<:Real}, B::Real = 0.0;
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





function write_cr_to_dat(t::Vector{<:Real}, CR::Vector{CRMomentumDistribution}, filename::String)

    Nsnaps = length(t)
    Nbins  = length(CR[1].norm)

    f = open(filename, "w")
    write(f, Nsnaps)
    write(f, Nbins)
    write(f, t)
    for i = 1:Nsnaps
        write(f, CR[i].bounds)
        write(f, CR[i].norm)
    end
    close(f)

end


function read_cr_from_dat(filename)

    f = open(filename, "r")
    Nsnaps = read(f, Int64)
    Nbins  = read(f, Int64)

    t = read!(f, Vector{Float64}(undef, Nsnaps))

    CR = Array{CRMomentumDistribution,1}(undef,Nsnaps)
    t = data[:,1]

    for i = 1:Nsnaps
        bounds = read!(f, Vector{Float64}(undef, Nbins+1))
        norm   = read!(f, Vector{Float64}(undef, Nbins))
        CR[i] = CRMomentumDistribution([data[i,2:2Nbins+1]; data[i,2Nbins+1]], data[i,2Nbins+2:4Nbins+1])
    end

    return t, CR
end

"""
    read_crp_cre_from_txt(filename::String)

Read CR Proton and Electron spectra from a txt file.
"""
function read_cr_from_txt(filename::String)

    data = readdlm(filename)

    N = size(data,1)
    Nbins = Int64((size(data,2) - 1) / 4)

    CR = Array{CRMomentumDistribution,1}(undef,N)

    t = data[:,1]

    for i = 1:N
        CR[i] = CRMomentumDistribution([data[i,2:2Nbins+1]; data[i,2Nbins+1]], data[i,2Nbins+2:4Nbins+1])
    end

    return t, CR
end


function write_cr_to_txt(t::Vector{<:Real}, CR::Vector{CRMomentumDistribution}, output_file::String)

    data = Matrix{Float64}(undef, length(t), 1+2*length(CR[1].norm))

    for i = 1:length(t)
        data[i,:] = [t[i] CR[i].bound[1:end-1]' CR[i].norm' ]
    end

    writedlm(output_file, data)
end
