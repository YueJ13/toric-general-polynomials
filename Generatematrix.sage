def lie_derivative_mono(monomial, g, n=0, d=0):
    '''
    Compute the Lie derivative of a monomial with respect to a given transformation matrix.
    Input:
    - monomial: A single-term polynomial whose Lie derivative is to be computed.
    - g: A matrix representing the infinitesimal transformation.
    - n: The number of variables in the ambient polynomial ring (optional, default=0).
    Output:
    - The resulting polynomial after applying the Lie derivative.
    '''
    if n == 0:
        n_vars = len(monomial.parent().gens())
    else:
        n_vars = n

    var = ['x%i' % i for i in range(1, n_vars + 1)]
    var += ['c%s' % ''.join(map(str, i)) for i in IntegerVectors(d, n_vars)]
    var += ['g%i%i' % (i, j) for i in range(1, n_vars + 1) for j in range(1, n_vars + 1)]
    B = PolynomialRing(QQ, var)
    B.inject_variables()

    monomial = B(monomial)
    
    X = matrix(B, n_vars, 1,  B.gens()[:n_vars])
    gen_index = {B.gens()[l]: l for l in range(len(B.gens()))}
    der = 0
    for x in monomial.factor():
        l = (g[gen_index[x[0]], :] * X)[0, 0]
        der += l * B(monomial / x[0]) * x[1]
    return der


def lie_derivative_poly(poly, g, n=0, d=0):
    '''
    Compute the Lie derivative of a polynomial with respect to a given transformation matrix.
    Input:
    - poly: The polynomial whose Lie derivative is to be computed.
    - g: A matrix representing the infinitesimal transformation.
    - n: The number of variables in the ambient polynomial ring (optional, default=0).
    Output:
    - The resulting polynomial after applying the Lie derivative.
    '''
    if n == 0:
        n_vars = len(poly.parent().gens())
    else:
        n_vars = n
    der = 0
    for mon, coef in zip(poly.monomials(), poly.coefficients()):
        der += lie_derivative_mono(mon, g, n_vars, d) * coef
    return der

def monomial_generator(variables, degree):
    '''
    Generate all monomials in the given variables and degree.
    Input:
    - variables: A list of variables.
    - degree: A positive integer representing the degree of the monomials.
    Output:
    - A list of all monomials in the given variables and degree, presented in the reverse lexicographic order.
    '''
    n_vars = len(variables)
    degs = list(WeightedIntegerVectors(degree, [1 for i in range(n_vars)]))
    monomials = []
    for d in degs:
        mon = 1
        for i in range(n_vars):
            mon *= variables[i]**d[i]
        monomials.append(mon)
    return monomials


def poly_to_vec(polynomial, n=0, degree=0):
    '''
    Convert a homogeneous polynomial into its vector representation based on a monomial basis.
    Input:
    - polynomial: The polynomial to be converted into a vector.
    - n: The number of variables in the ambient polynomial ring (optional, default=0).
    - degree: The degree of the polynomial (optional, default=0).
    Output:
    - A list representing the vectorized form of the polynomial, where each entry corresponds to the coefficient of a monomial in the chosen basis.
    '''
    if n == 0:
        n_vars = len(R.gens())
    else:
        n_vars = n
    if degree == 0:
        d = polynomial.degree()
    else:
        d = degree

    var = ['x%i' % i for i in range(1, n_vars + 1)]
    var += ['c%s' % ''.join(map(str, i)) for i in IntegerVectors(d, n_vars)]
    var += ['g%i%i' % (i, j) for i in range(1, n_vars + 1) for j in range(1, n_vars + 1)]
    R = PolynomialRing(QQ, var)
    R.inject_variables()

    n_2 = binomial(n_vars + d - 1, d)
    monomials = monomial_generator(R.gens()[:n_vars], d)
    vec = [0] * n_2
    for i in range(n_2):
        vec[i] = polynomial.coefficient(monomials[i])
    return vec

def symmalg_homo(generators, g, n=0, d=0):
    '''
    Compute the "A" matrix of homogeneous polynomials.
    Input:
    - generators: A list of generating polynomials defining the homogeneous ideal. All generators must be homogeneous polynomials of the same degree.
    - n: The number of variables in the ambient polynomial ring (optional, default=0).
    Output:
    - A: The coefficient matrix of (g11, ..., gnn, lambda_1, ..., lambda_k), where k is the number of generators. 
    '''
    S = generators[0].parent()
    if n != 0:
        n_vars = n
    else:
        n_vars = len(S.gens())
    
    der_gens_vecs = [poly_to_vec(lie_derivative_poly(gen, g, n_vars, d), n_vars, d) for gen in generators]
    # print('\nvector of lie derivative:\n')
    # print(der_gens_vecs)
    
    var = ['x%i' % i for i in range(1, n_vars + 1)]
    var += ['c%s' % ''.join(map(str, i)) for i in IntegerVectors(d, n_vars)]
    var += ['g%i%i' % (i, j) for i in range(1, n_vars + 1) for j in range(1, n_vars + 1)]
    R = PolynomialRing(QQ, var)
    R.inject_variables()
    gens_vecs = [poly_to_vec(R(str(gen)), n_vars, d) for gen in generators]
    # print('\nvector of generators:\n')
    # print(gens_vecs)

    ker = []
    for eqn in der_gens_vecs[0]:
        l = [eqn.coefficient(gij) for gij in R.gens()[n_vars+len(IntegerVectors(d, n_vars)):]]
        ker.append(l)
    K = matrix(ker)

    A = K.augment(vector(gens_vecs[0]).column())
    print('\nmatrix A:\n')
    print(A)