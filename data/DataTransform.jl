"""
######################## NETWORK TRANSFORM #########################
Functions to turn a single phase file into a three phases.
"""
module DataTransform

function from_single_phase_to_three_phases(inpath::String,outpath::String)
    lines = []  # Create an array to store lines
    open(inpath, "r") do file
        for line in eachline(file)
            push!(lines, strip(line))  # Add each line to the array after stripping whitespace
        end
    end
    # Procedure to transform the lines
    open(outpath, "w") do file
        for (i,line) in enumerate(lines)
            # Change phases from 1 to 3
            if occursin("hases", line)
                line = replace(line, r"hases=\d+(\.\d+)?" => "hases=3")
            end
            # Change buses from 1 to 3
            if occursin("Bus", line)
                line = replace(line, r"(Bus\d+\=\d+\.\d+)" => s"\1.2.3")
            end
            # Change linecodes from 1 to 3
            if occursin("linecode", line)
                line = replace(line, r"(nphases=\d+)" => s"\1 basefreq=50")
                # Find impedance values
                r_value_match = match(r"Rmatrix=\[(\d+\.\d+)\]", line); r_value = parse(Float64, r_value_match[1])
                x_value_match = match(r"Xmatrix=\[(\d+\.\d+)\]", line); x_value = parse(Float64, r_value_match[1])
                # Define impedance values
                r_matrix = [r_value (0.04*0.1)/r_value r_value (0.04*0.1)/r_value (0.04*0.1)/r_value r_value]; r_matrix_str = "Rmatrix=[" * format_pairs(r_matrix) * "]"
                x_matrix = [x_value (0.02*0.1)/x_value x_value (0.02*0.1)/x_value (0.02*0.1)/x_value x_value]; x_matrix_str = "\n ~ Xmatrix=[" * format_pairs(x_matrix) * "]"
                c_matrix = [51 -0 51 -0 -0 51]; c_matrix_str = "\n ~ Cmatrix=[" * format_pairs(c_matrix) * "]"
                # Place impedance values
                line = replace(line, r"(Rmatrix=\[(\d+\.\d+)\])" => r_matrix_str)
                line = replace(line, r"Xmatrix=\[(\d+\.\d+)\]" => x_matrix_str)
                line *= c_matrix_str
            end
            write(file, line * "\n")  # Write each line with a newline
        end
    end
end

function from_three_phases_to_single_phase(inpath::String,outpath::String)
    lines = []  # Create an array to store lines
    open(inpath, "r") do file
        for line in eachline(file)
            push!(lines, strip(line))  # Add each line to the array after stripping whitespace
        end
    end
    # Procedure to transform the lines
    open(outpath, "w") do file
        for (i,line) in enumerate(lines)
            # Change phases from 1 to 3
            if occursin("hases", line)
                line = replace(line, r"hases=\d+(\.\d+)?" => "hases=1")
            end
            # Change buses from 1 to 3
            if occursin("Bus", line)
                line = replace(line, r"(Bus\d+\=\d+\.\d+)" => s"\1")
            end
            # Change linecodes from 1 to 3
            if occursin("linecode", line)
                line = replace(line, r"(nphases=\d+)" => s"\1 basefreq=50")
                # Find impedance values
                r_value_match = match(r"Rmatrix=\[(\d+\.\d+)\]", line); r_value = parse(Float64, r_value_match[1])
                x_value_match = match(r"Xmatrix=\[(\d+\.\d+)\]", line); x_value = parse(Float64, r_value_match[1])
                # Define impedance values
                r_matrix = [r_value]; r_matrix_str = "Rmatrix=[" * r_matrix * "]"
                x_matrix = [x_value]; x_matrix_str = "\n ~ Xmatrix=[" * x_matrix * "]"
                c_matrix = [51]; c_matrix_str = "\n ~ Cmatrix=[" * c_matrix * "]"
                # Place impedance values
                line = replace(line, r"(Rmatrix=\[(\d+\.\d+)\])" => r_matrix_str)
                line = replace(line, r"Xmatrix=\[(\d+\.\d+)\]" => x_matrix_str)
                line *= c_matrix_str
            end
            write(file, line * "\n")  # Write each line with a newline
        end
    end
end

end