# mtimesx

> Fast Matrix Multiply with Multi-Dimensional Support

**MTIMESX** is a fast general purpose matrix and scalar multiply routine that has the following features:

- Supports multi-dimensional (nD, n>2) arrays directly
- Supports Transpose, Conjugate Transpose, and Conjugate pre-operations
- Supports singleton expansion
- Utilizes BLAS calls, custom C loop code, or OpenMP multi-threaded C loop code
- Can match MATLAB results exactly or approximately as desired
- Can meet or beat MATLAB for speed in most cases

**MTIMESX** has six basic operating modes:

- **BLAS** - Always uses BLAS library calls
- **LOOPS** - Always uses C loops if available
- **LOOPSOMP** - Always uses OpenMP multi-threaded C loops if available
- **MATLAB** -, Fastest BLAS or LOOPS method that matches MATLAB exactly (default)
- **SPEED** - Fastest BLAS or LOOPS method even if it doesn't match MATLAB exactly
- **SPEEDOMP** - Fastest BLAS, LOOPS, or LOOPOMP method even if it doesn't match MATLAB exactly

**MTIMESX** inputs can be:

- single
- double
- double sparse

## Getting Started

The general syntax is (arguments in brackets `[]` are optional):

```
mtimesx( [directive] )
mtimesx( A [,transa] ,B [,transb] [,directive] )
```

Where `transa`, `transb`, and `directive` are the optional inputs:

- `transa` - A character indicating a pre-operation on A
- `transb` - A character indicating a pre-operation on B

    The pre-operation can be any of:

    - `N` or `n` - No pre-operation (the default if trans is missing) 
    - `T` or `t` - Transpose
    - `C` or `c` - Conjugate Transpose
    - `G` or `g` - Conjugate (no transpose)

- `directive` - One of the operation modes(.i.e. **BLAS** or **LOOPS**) listed above, or other directives

### Examples:

```
C = mtimesx(A,B) % performs the calculation C = A * B
C = mtimesx(A,'T',B) % performs the calculation C = A.' * B
C = mtimesx(A,B,'g') % performs the calculation C = A * conj(B)
C = mtimesx(A,'c',B,'C') % performs the calculation C = A' * B'
mtimesx('SPEEDOMP','OMP_SET_NUM_THREADS(4)') % sets SPEEDOMP mode with number of threads = 4
```

For nD cases, the first two dimensions specify the matrix multiply involved. The remaining dimensions are duplicated and specify the number of individual matrix multiplies to perform for the result. i.e., MTIMESX treats these cases as arrays of 2D matrices and performs the operation on the associated parings. For example:

If `A` is `(2,3,4,5)` and `B` is `(3,6,4,5)`, then `mtimesx(A,B)` would result in `C(2,6,4,5)`, where `C(:,:,i,j) = A(:,:,i,j) * B(:,:,i,j)`, `i=1:4`, `j=1:5`

which would be equivalent to the MATLAB m-code:

```matlab
C = zeros(2,6,4,5);
for m=1:4
    for n=1:5
            C(:,:,m,n) = A(:,:,m,n) * B(:,:,m,n);
    end
end
```

The first two dimensions must conform using the standard matrix multiply rules taking the transa and transb pre-operations into account, and dimensions 3:end must match exactly or be singleton (equal to 1). If a dimension is singleton then it is virtually expanded to the required size (i.e., equivalent to a repmat operation to get it to a conforming size but without the actual data copy). This is equivalent to a bsxfun capability for matrix multiplication.
