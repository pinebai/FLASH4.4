function add_var, input_names, new_name

new_size = (size(input_names))[1] + 1
new_varnames = strarr(new_size)


new_varnames[0:new_size - 2] = input_names
new_varnames[new_size - 1]   = new_name

return, new_varnames
end
