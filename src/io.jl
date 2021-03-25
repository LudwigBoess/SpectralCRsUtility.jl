using GadgetUnits
using GadgetIO

global const cnst_c = 2.9979e10
global const slope_soft = 1.e-6

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


function construct_spectrum(CR_N, CR_S, CR_C, pmin, pmax, mc)
    
    Nbins = size(CR_N,1)
    bin_width = log10(pmax/pmin) / Nbins

    CR_dis   = Array{Float64,1}(undef, Int(2*Nbins))
    CR_bound = Array{Float64,1}(undef, Int(2*Nbins + 1))

    # get zeroth bin
    CR_bound[1] = pmin
    CR_dis[1] = CR_N[1]

    # all other bins
    j = 2
    for i = 1:Nbins-1

        # upper boundary of bin
        CR_bound[j] = pmin * 10.0^(bin_width*i)
        CR_dis[j] = CR_N[i] * ( CR_bound[j]/CR_bound[j-1])^(-CR_S[i])
        if CR_bound[j] > CR_C/mc
            CR_bound[j] = CR_C/mc
        end

        # lower bound of next bin
        CR_bound[j+1] = CR_bound[j]
        CR_dis[j+1] = CR_N[i+1]
        if CR_bound[j] == CR_C/mc
            CR_bound[j+1] = CR_C/mc
        end

        j += 2

    end

    # last boundary
    CR_bound[j] = pmax
    if CR_bound[j-1] < CR_C/mc
        CR_bound[j] = CR_C/mc
    end
    CR_dis[j]     = CR_N[Nbins] * ( CR_bound[j]/CR_bound[j-1])^(-CR_S[Nbins])
    CR_bound[j+1] = CR_bound[j]

    return CR_bound, CR_dis

end

function getCRMomentumDistributionFromPartID(snap_file::String, ID::Integer;
                                             pmin::Real=1.0, pmax::Real=1.0e6,
                                             Nbins::Integer=0, mode::Int64=2)

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

    # electrons
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

    Nbins = size(CRpS,1)

    par = CRMomentumDistributionConfig(pmin, pmax, Nbins, mode)

    cr = CRMomentumDistribution(par.Nbins)

    # compensate for io - needs to be converted to Float64!
    if mode == 1
        CRpN = 1.e20 .* Float64.(CRpN)
        CReN = 1.e20 .* Float64.(CReN)
    elseif mode == 2
        CRpN = 1.e-20 .* Float64.(CRpN)
        CReN = 1.e-20 .* Float64.(CReN)
    elseif mode == 3
        @inbounds for i = 1:Nbins
            CRpN[i] = 10.0^CRpN[i]
            CReN[i] = 10.0^CReN[i]
        end
    end

    cr.CRp_bound, cr.CRp_dis = construct_spectrum(CRpN, CRpS, CRpC, pmin, pmax, par.mc_p)
    cr.CRe_bound, cr.CRe_dis = construct_spectrum(CReN, CReS, CReC, pmin, pmax, par.mc_e)

    return cr
end