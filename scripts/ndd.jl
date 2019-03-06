#! /usr/local/julia-1.1.0/bin/julia
###############################################################################
# Utility to displace a PDB along many vectors, generating 1 PDB for each
# vector, useful for NDD analysis.
# Calpha only.
# by https://github.com/pgbarletta
###############################################################################
using Chemfiles, JUMD
using ArgParse
using DelimitedFiles, LinearAlgebra, FileIO

function displaceAA(in_frm, aa, aa_3, in_vec)
    # Preparo variables
    in_top = Topology(in_frm)
    natoms = convert(Int64, size(in_top))
    in_xyz = positions(in_frm)

    # Determino orden de residuos (hay q actualizar el Julia Chemfiles)
    tmp = Array{Int64}(undef, aa)
    ids = Array{Int64}(undef, aa)
    [ ids[i+1] = convert(Int64, id((Residue(in_top, i)))) for i = 0:aa-1 ]
    idx = sortperm(ids)
    # Determino el nro de atomos de c/ aminoácido
    [ tmp[i+1] = size(Residue(in_top, i)) for i = 0:aa-1 ]
    natom_aa = tmp[idx]

    # Paso el vector columna de tamaño 1xaa_3 a 3xaa
    vector = reshape(in_vec, 3, aa)
    # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    sum_mat = Array{Float64}(undef, 3, natoms)
    cursor = 0
    for i = 1:aa
        if i == 1
            sum_mat[:, 1:natom_aa[i]] = repeat(vector[:, 1], 1, natom_aa[i])
            cursor = natom_aa[i]
            continue
        end
        rango = collect(cursor+1:cursor + natom_aa[i])
        sum_mat[:, rango] = repeat(vector[:, i], 1, natom_aa[i])
        cursor += natom_aa[i]
    end

    # Listo, ahora puedo mover el pdb
    out_frm = deepcopy(in_frm)
    out_xyz = positions(out_frm)

    # Tengo q hacer esto por ahora, hasta q arreglemos Chemfiles.
    for i = 1:size(in_xyz)[1]
        for j = 1:size(in_xyz)[2]
            out_xyz[i, j]  = round(in_xyz[i, j] + sum_mat[i, j], digits = 3)
        end
    end
    return out_frm
end

# Arg Parse settings
s = ArgParseSettings()
@add_arg_table s begin
    "--in_pdb_filename", "-p"
        help = "Input PDB."
        arg_type = String
        required = true
    "--modes_filename", "-v"
        help = "Input modes."
        arg_type = String
        required = true
    "--mul", "-m"
        help = "Multiplier."
        arg_type = Int
        required = true
    "--suffix", "-o"
        help = "Output PDBs suffix"
        arg_type = String
        required = true
    "--out_dir", "-d"
        help = "Output PDBs directory"
        arg_type = String
        required = true
    "--amber_modes", "-a"
        help = "Mark true when reading Amber(format) modes. Default: false."
        arg_type = Bool
        required = false
        default = false
    "--weights_filename", "-w"
        help = "Input 1/weights, if desired. If Amber modes are read, eigenvalues will be used. Default: none"
        arg_type = String
        required = false
        default = "none"
end

# Read arguments from console
parsed_args = parse_args(ARGS, s)
args = Array{Any, 1}(undef, 0)
for (arg, val) in parsed_args
    arg = Symbol(arg)
    @eval (($arg) = ($val))
end

println("Input parameters:")
println("INPDB          ", in_pdb_filename)
println("MODES          ", modes_filename)
println("MUL            ", mul)
println("SUFFIX         ", suffix)
println("OUT_DIR        ", out_dir)
println("AMBER_MODES    ", amber_modes)
println("WEIGHTS        ", weights_filename)

# Read PDB
const in_trj = Trajectory(in_pdb_filename)
const in_frm = read(in_trj)
const in_top = Topology(in_frm)
const aa = convert(Int64, count_residues(in_top))
const aa_3 = aa * 3
nmodes = aa_3 - 6

# Read input vectors and weights
in_modes = Array{Float64}(undef, aa_3, nmodes)
weights = Array{Float64}(undef, nmodes)

if (amber_modes)
    try
        global in_modes, weights = JUMD.readPtrajModes(modes_filename)
    catch e
        println("Error when reading Amber modes: ", modes_filename)
        error(e)
    end
else
    try
        global in_modes = convert(Array{Float64}(aa_3, nmodes), readdlm(modes_filename))
    catch e
        println("Error when reading modes: ", modes_filename)
        error(e)
    end
    if weights_filename != "none"
        try
            global weights = convert(Array{Float64}(nmodes), readdlm(weights_filename))
        catch e
            println("Error when reading weights: ", weights_filename)
            error(e)
        end
    else
        global weights = convert(Array{Float64}(nmodes), fill(1., nmodes))
    end
end


# Ahora desplazo
pdb_names = Array{String}(undef, nmodes)
for i = 1:nmodes
    # Escalo vector
    modo = in_modes[:, i] .* mul ./ weights[i]
    # Desplazo
    out_frm = displaceAA(in_frm, aa, aa_3, modo);
    # Y guardo
    pdb_names[i] = joinpath(out_dir, string(i, "_", suffix, ".pdb"))
    out_trj = Trajectory(pdb_names[i], 'w')
    write(out_trj, out_frm)
    close(out_trj)
end

# Write in_ndd file
const in_ndd_filename = string("in_ndd_",
    splitext(basename(in_pdb_filename))[1])
writedlm(joinpath(out_dir, in_ndd_filename), pdb_names)
