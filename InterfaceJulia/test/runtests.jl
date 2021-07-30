# Empty file for testing Julia code
using InterfaceJulia.TasData

# Test to create a lower set of indices of total degree 3
foo = create_lower_set(3, x -> TasData.is_itd_elem(5, x))
