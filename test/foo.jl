using ArgParse

argset = ArgParseSettings()
@add_arg_table argset begin
    "-r"
    "-s"
end
args = parse_args(argset)

r = args["r"]
s = args["s"]

print("$r/$s\n")