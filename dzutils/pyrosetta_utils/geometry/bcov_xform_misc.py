def xform_to_rotation(xform):
    mat = numeric.xyzMatrix_double_t()
    mat.xx = xform[_xx]
    mat.xy = xform[_xy]
    mat.xz = xform[_xz]
    mat.yx = xform[_yx]
    mat.yy = xform[_yy]
    mat.yz = xform[_yz]
    mat.zx = xform[_zx]
    mat.zy = xform[_zy]
    mat.zz = xform[_zz]

    return mat


def rotation_to_xform(mat):
    xform = np.identity(4)
    xform[_xx] = mat.xx
    xform[_xy] = mat.xy
    xform[_xz] = mat.xz
    xform[_yx] = mat.yx
    xform[_yy] = mat.yy
    xform[_yz] = mat.yz
    xform[_zx] = mat.zx
    xform[_zy] = mat.zy
    xform[_zz] = mat.zz

    return xform


def xform_to_vector(xform):
    vec = numeric.xyzVector_double_t()
    vec.x = xform[_x]
    vec.y = xform[_y]
    vec.z = xform[_z]

    return vec


def vector_to_xform(vec):
    xform = np.identity(4)
    xform[_x] = vec.x
    xform[_y] = vec.y
    xform[_z] = vec.z

    return vec
