# BSP
A c++ bsp-tree for 3d triangle meshes that can handle arbitrary vertex types. Look at the examples for use-cases and at the hpp file for documentation

the 3 examples are all relatively similar but I tried to include some variation

I did try to make it as flexible as possible regarding the vertex type. The followind constraints apply:

- one vertex must be contained within a datatype (struct, array of bytes, string containing descritpion, index into another array)
- that datatype must be possible to put into a vector
- you must provide a getter-method to get the 3d-position of that vertex (using template specialization)
- the returntype of that function must be boost:qvm compatible (which can be acieved using template specialization)
- you must provide a function (via template specialization) that adds a new vertex at the end of a vector
