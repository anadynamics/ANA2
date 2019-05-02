using DelimitedFiles, LinearAlgebra, FileIO, DataFrames, CSV

function get_opciones(file::String)
# Isolate lines corresponding to options.
    sz = length(file)
    lines = Array{String, 1}()
    let i = 1
        while i < sz - 5
            if file[i:i+5] == "    (\""
                let j = i + 6
                    for j = i+6:1:sz-10
                        if file[j:j+3] == "\\n\")"
                            push!(lines, file[i+6:j-2])
                            i = j + 2
                            break
                        end
                    end
                end
            end
            i += 1
        end
    end

    return lines
end

function get_keyword(linea::AbstractString)
    var = ""
    def_val = ""
    let var, def_val, a = match(r"default_value\(", linea).offset
        for i = 1:length(linea)
            if linea[i] == '&'
                for j = i:length(linea)
                    if linea[j] == ')'
                        var = linea[i+1:j-1]
                        break
                    end
                end
                break
            end
        end
        if a != Nothing
            a += 14
            b = findfirst(isequal(')'), linea[a:end])
            def_val = linea[a:a+b-2]
        end
        return var, replace(def_val, "\"" => "")
    end
end

function padea(v::Array{String, 1}, n::Int64)
    n = length(v)
    while n < tope
        push!(v, "")
        n+=1
    end
end

#########################
# Script medio chorro p/ sacar la documentación de ANA del ProgramOptions.cpp
# No hay q poner comas (",") en el texto de config pq sino se caga todo.
#########################

io = open("../src/ProgramOptions.cpp", "r")
file = read(io, String)
close(io)
raw_opciones = get_opciones(file)



opciones = Array{String, 1}()
ops = Array{String, 1}()
textos = Array{String, 1}()
default_vars = Array{String, 1}()
variables = Array{String, 1}()

let terminal_opts = true
    for i = 1:length(raw_opciones)
        tempo = split(raw_opciones[i], ",")
        # Remove weird chars
        texto = replace(replace(replace(tempo[end], "\"" => ""), "\\n" => ""),
        "\n" => "")
        push!(textos, texto)
        linea_var = String
        if terminal_opts
            push!(opciones, tempo[1])
            push!(ops, string(tempo[2][1]))
            global linea_var = tempo[3]
        else
            push!(opciones, tempo[1][1:end-1])
            push!(ops, "-")
            global linea_var = tempo[2]
        end

        if opciones[i] == "help"
            terminal_opts = false
        end
        if opciones[i] == "help" || opciones[i] == "ver"
            push!(variables, "-")
            push!(default_vars, "-")
            continue
        end
        var, def_val = get_keyword(linea_var)
        push!(variables, var)
        push!(default_vars, def_val)
    end
end

df_opts = DataFrame(Name = opciones, Shorthand = ops, Default = default_vars,
    Description = textos, Variables = variables)
CSV.write("help.csv", df_opts, delim = '\t')
