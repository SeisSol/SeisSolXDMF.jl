module SeisSolXDMF

export XDMFFile, timesteps_of, grid_of, data_of

using XMLDict
using HDF5

    struct XDMFFile
        xml         :: XMLDict.XMLDictElement
        base_path   :: AbstractString
        timesteps   :: Array{XMLDict.XMLDictElement,1}
    end

    function expect(expression, error_msg)
        if !expression
            error(error_msg)
        end
    end

    function XDMFFile(file_path::AbstractString)
        xml_file = open(file_path) do file; read(file, String) end
        xml_file = parse_xml(xml_file)

        try
            # The root node "Xdmf" is assumed
            timesteps = xml_file["Domain"]["Grid"]["Grid"]
            num_timesteps = length(timesteps)

            expect(num_timesteps >  0, "No timesteps found in XDMF file.")

            xdmf_base_path = splitdir(file_path)[1]
            return XDMFFile(xml_file, xdmf_base_path, timesteps)
        catch
            error("Input file has no valid SeisSol XDMF structure.")
        end
    end

    function timesteps_of(xdmf::XDMFFile) :: Array{Float64, 1}
        return Array(map(data_item -> parse(Float64, data_item["Time"][:Value]), xdmf.timesteps))
    end

    function grid_of(xdmf::XDMFFile) :: Tuple{AbstractArray{Integer, 2}, AbstractArray{AbstractFloat, 2}}
        # Read geometry (points) and topology (simplices (triangles/tetrahedra)) files.
        # Since the ids of the points referred to in the topology array start at 0 (and Julia counts from 1) we have to add 1.
        simplices = read_dataset(xdmf, xdmf.timesteps[1]["Topology"]["DataItem"]) .+ 1
        expect(ndims(simplices) == 2, "Topology should have 2 dimensions but has $(ndims(simplices)).")
        expect(size(simplices, 1) ∈ (3, 4), "Topology entries must have 3 (triangles) or 4 (tetrahedra) points but have $(size(simplices, 1)).")

        points = read_dataset(xdmf, xdmf.timesteps[1]["Geometry"]["DataItem"])
        expect(ndims(points) == 2, "Geometry should have 2 dimensions but has $(ndims(simplices)).")
        expect(size(points, 1) == 3, "Geometry entries must have 3 coordinates but have $(size(simplices, 1)).")

        return (simplices, points)
    end

    function data_of(xdmf::XDMFFile, timestep::Integer, var_name::AbstractString)
        expect(timestep >= 1, "Timestep $(timestep) is smaller than 1.")
        expect(timestep <= length(xdmf.timesteps), "Timestep $(timestep) exceeds number of timesteps in XDMF file ($(length(xdmf.timesteps))).")

        attrs = xdmf.timesteps[timestep]["Attribute"]
        var_index = findfirst(attr -> attr[:Name] == var_name, attrs)
        expect(!isnothing(var_index), """Variable $var_name not found in timestep $timestep. 
            Found variables: $(join(map(attr -> attr[:Name], attrs), ", ")).""")

        var_entry = attrs[var_index]["DataItem"]

        return reshape(read_hyperslab(xdmf, var_entry), :)
    end
        
    function read_dataset(xdmf::XDMFFile, data_item::XMLDict.XMLDictElement; indices=nothing)
        filename = data_item[""]

        file_range = data_item[:Dimensions]
        file_range = split(file_range, ' ')
        file_range = parse.(UInt, file_range)
        file_period = file_range[2]

        number_type = get_number_type(data_item)

        if data_item[:Format] == "HDF"
            path_parts = split(filename, ':', limit=2)
            expect(!isempty(path_parts[1]) && !isempty(path_parts[2]), "HDF5 group path is invalid.")

            filename = joinpath(xdmf.base_path, path_parts[1])
            hdf_dataset_path = String(path_parts[2])

            dataset = 
                if isnothing(indices)
                    h5read(filename, hdf_dataset_path)
                else
                    h5read(filename, hdf_dataset_path, (
                        indices[1,2]:indices[2,2]:(indices[1,2] + indices[3,2] - 1),
                        indices[1,1]:indices[2,1]:(indices[1,1] + indices[3,1] - 1)))
                end
        else
            filename = joinpath(xdmf.base_path, filename)
            start_index = (indices[1,1]-1) * file_period + (indices[1,2]-1) # Starting at 0
            stride = (indices[2,1]-1) * file_period + (indices[2,2]-1)
            expect(stride == 0, "Stride not yet supported for POSIX output.") # (no need (?))
            dataset = Array{number_type, 2}(undef, (indices[3,1], indices[3,2]))
            open(filename, "r") do file
                seek(file, start_index * sizeof(number_type))
                read!(file, dataset)
            end
        end

        return dataset
    end

    function read_hyperslab(xdmf::XDMFFile, data_item::XMLDict.XMLDictElement)
        #=
        Example DataItem:
        <DataItem ItemType="HyperSlab" Dimensions="4561037">
            <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 0 1 1 1 4561037</DataItem>
            <DataItem NumberType="Float" Precision="8" Format="Binary" Dimensions="1 5242880">out_cell/mesh0/u.bin</DataItem>
        </DataItem>
        =#
        expect(data_item[:ItemType] == "HyperSlab", "Unexpected item type '$(data_item[:ItemType])', expected 'HyperSlab'.")
        range_item = data_item["DataItem"][1]
        binary_item = data_item["DataItem"][2]

        expect(range_item[:Format] == "XML", "Unexpected format '$(range_item[:Format])', expected 'XML'.")
        expect(range_item[:Dimensions] == "3 2", "Unexpected dimensions '$(range_item[:Dimensions])', expected '3 2'.")
        expect(binary_item[:Format] ∈ ["Binary", "HDF"], "Unexpected format '$(binary_item[:Format])', expected 'Binary' or 'HDF'.")

        rng = range_item[""]
        rng = split(rng, ' ')
        rng = parse.(Int, rng)
        rng = reshape(rng, (2, 3))'

        rng[1,:] .+= 1 # Julia start indices

        return read_dataset(xdmf, binary_item, indices=rng)
    end

    function get_number_type(data_item::XMLDict.XMLDictElement) :: Type
        bytes_of_precision = parse(UInt, data_item[:Precision])
        number_type_name = data_item[:NumberType]

        print_error(type_name, precision) = 
            error("""Type "$type_name" of size $precision bytes is not recognized by the parser.""")

        if bytes_of_precision == 8
            if number_type_name == "Int"
                return Int64
            elseif number_type_name == "UInt"
                return UInt64
            elseif number_type_name == "Float"
                return Float64
            else
                print_error(number_type_name, bytes_of_precision)
            end
        elseif bytes_of_precision == 4
            if number_type_name == "Int"
                return Int32
            elseif number_type_name == "UInt"
                return UInt32
            elseif number_type_name == "Float"
                return Float32
            else
                print_error(number_type_name, bytes_of_precision)
            end
        else
            print_error(number_type_name, bytes_of_precision)
        end
    end

end # module
