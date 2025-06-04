function [x, w] = GJWT ( order, alf, bt, varargin )
%*****************************************************************************80
%
%% jacobi_rule() generates a Gauss-Jacobi rule.
%
%  Discussion:
%
%    This program computes a standard Gauss-Jacobi quadrature rule
%    and writes it to a file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    28 February 2010
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    ORDER (number of points) in the rule;
%
%    ALPHA, the exponent of (1+x);
%
%    BETA, the exponent of (1-x);
%    
%    A, the left endpoint;
%
%    B, the right endpoint;
%
%    FILENAME, the root name of the output files.
%
    if length(varargin)==2
        a = varargin{1};
        b = varargin{2};
    else
        a = 0;
        b = 1;
    end
    kind = 4;
    [ x, w ] = cgqf ( order, kind, bt-1, alf-1, a, b );
end

function [ t, wts ] = cdgqf ( nt, kind, bt, alf )

%*****************************************************************************80
%
%% cdgqf() computes a Gauss quadrature formula with default A, B and simple knots.
%
%  Discussion:
%
%    This routine computes all the knots and weights of a Gauss quadrature
%    formula with a classical weight function with default values for A and B,
%    and only simple knots.
%
%    There are no moments checks and no printing is done.
%
%    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 January 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Input:
%
%    integer NT, the number of knots.
%
%    integer KIND, the rule.
%    1, Legendre,             (a,b)       1.0
%    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
%    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^bt
%    4, Jacobi,               (a,b)       (b-x)^bt*(x-a)^alf
%    5, Generalized Laguerre, (a,inf)     (x-a)^bt*exp(-b*(x-a))
%    6, Generalized Hermite,  (-inf,inf)  |x-a|^bt*exp(-b*(x-a)^2)
%    7, Exponential,          (a,b)       |x-(a+b)/2.0|^bt
%    8, Rational,             (a,inf)     (x-a)^bt*(x+b)^alf
%
%    real BETA, the value of Alpha, if needed.
%
%    real ALPHA, the value of Beta, if needed.
%
%  Output:
%
%    real T(NT), the knots.
%
%    real WTS(NT), the weights.
%
  parchk ( kind, 2 * nt, bt, alf );
%
%  Get the Jacobi matrix and zero-th moment.
%
  [ aj, bj, zemu ] = class_matrix ( kind, nt, bt, alf );
%
%  Compute the knots and weights.
%
  [ t, wts ] = sgqf ( nt, aj, bj, zemu );

  return
end
function [ t, wts ] = cgqf ( nt, kind, bt, alf, a, b )

%*****************************************************************************80
%
%% cgqf() computes knots and weights of a Gauss quadrature formula.
%
%  Discussion:
%
%    The user may specify the interval (A,B).
%
%    Only simple knots are produced.
%
%    The user may request that the routine print the knots and weights,
%    and perform a moment check.
%
%    Use routine EIQFS to evaluate this quadrature formula.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 February 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Input:
%
%    integer NT, the number of knots.
%
%    integer KIND, the rule.
%    1, Legendre,             (a,b)       1.0
%    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
%    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^bt
%    4, Jacobi,               (a,b)       (b-x)^bt*(x-a)^alf
%    5, Generalized Laguerre, (a,+oo)     (x-a)^bt*exp(-b*(x-a))
%    6, Generalized Hermite,  (-oo,+oo)   |x-a|^bt*exp(-b*(x-a)^2)
%    7, Exponential,          (a,b)       |x-(a+b)/2.0|^bt
%    8, Rational,             (a,+oo)     (x-a)^bt*(x+b)^alf
%    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
%
%    real BETA, the value of Alpha, if needed.
%
%    real ALPHA, the value of Beta, if needed.
%
%    real A, B, the interval endpoints.
%
%  Output:
%
%    real T(NT), the knots.
%
%    real WTS(NT), the weights.
%

%
%  Compute the Gauss quadrature formula for default values of A and B.
%
  [ t, wts ] = cdgqf ( nt, kind, bt, alf );
%
%  All knots have multiplicity = 1.
%
  mlt = zeros(nt,1);
  mlt(1:nt) = 1;
%
%  NDX(I) = I.
%
  ndx = ( 1 : nt );
%
%  Scale the quadrature rule.
%
  [ t, wts ] = scqf ( nt, t, mlt, wts, nt, ndx, kind, bt, alf, a, b );

  return
end
function [ aj, bj, zemu ] = class_matrix ( kind, m, bt, alf )

%*****************************************************************************80
%
%% class_matrix() computes the Jacobi matrix for a quadrature rule.
%
%  Discussion:
%
%    This routine computes the diagonal AJ and subdiagonal BJ
%    elements of the order M tridiagonal symmetric Jacobi matrix
%    associated with the polynomials orthogonal with respect to
%    the weight function specified by KIND.
%
%    For weight functions 1-7, M elements are defined in BJ even
%    though only M-1 are needed.  For weight function 8, BJ(M) is
%    set to zero.
%
%    The zero-th moment of the weight function is returned in ZEMU.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 January 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Input:
%
%    integer KIND, the rule.
%    1, Legendre,             (a,b)       1.0
%    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
%    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^bt
%    4, Jacobi,               (a,b)       (b-x)^bt*(x-a)^alf
%    5, Generalized Laguerre, (a,inf)     (x-a)^bt*exp(-b*(x-a))
%    6, Generalized Hermite,  (-inf,inf)  |x-a|^bt*exp(-b*(x-a)^2)
%    7, Exponential,          (a,b)       |x-(a+b)/2.0|^bt
%    8, Rational,             (a,inf)     (x-a)^bt*(x+b)^alf
%
%    integer M, the order of the Jacobi matrix.
%
%    real BETA, the value of Alpha, if needed.
%
%    real ALPHA, the value of Beta, if needed.
%
%  Output:
%
%    real AJ(M), BJ(M), the diagonal and subdiagonal
%    of the Jacobi matrix.
%
%    real ZEMU, the zero-th moment.
%
  temp = eps;

  parchk ( kind, 2 * m - 1, bt, alf );

  temp2 = 0.5;

  if ( 500.0 * temp < abs ( ( gamma ( temp2 ) )^2 - pi ) )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'CLASS - Fatal error!\n' );
    fprintf ( 1, '  Gamma function does not match machine parameters.\n' );
    error ( 'CLASS - Fatal error!' );
  end

  bj = zeros(m,1);
  aj = zeros(m,1);

  if ( kind == 1 )

    ab = 0.0;

    zemu = 2.0 / ( ab + 1.0 );

    aj(1:m) = 0.0;

    for i = 1 : m
      abi = i + ab * mod ( i, 2 );
      abj = 2 * i + ab;
      bj(i) = abi * abi / ( abj * abj - 1.0 );
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 2 )

    zemu = pi;

    aj(1:m) = 0.0;

    bj(1) =  sqrt ( 0.5 );
    bj(2:m) = 0.5;

  elseif ( kind == 3 )

    ab = bt * 2.0;
    zemu = 2.0^( ab + 1.0 ) * gamma ( bt + 1.0 )^2 ...
      / gamma ( ab + 2.0 );

    aj(1:m) = 0.0;
    bj(1) = 1.0 / ( 2.0 * bt + 3.0 );
    for i = 2 : m
      bj(i) = i * ( i + ab ) / ( 4.0 * ( i + bt )^2 - 1.0 );
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 4 )

    ab = bt + alf;
    abi = 2.0 + ab;
    zemu = 2.0^( ab + 1.0 ) * gamma ( bt + 1.0 ) ...
      * gamma ( alf + 1.0 ) / gamma ( abi );
    aj(1) = ( alf - bt ) / abi;
    bj(1) = 4.0 * ( 1.0 + bt ) * ( 1.0 + alf ) ...
      / ( ( abi + 1.0 ) * abi * abi );
    a2b2 = alf * alf - bt * bt;

    for i = 2 : m
      abi = 2.0 * i + ab;
      aj(i) = a2b2 / ( ( abi - 2.0 ) * abi );
      abi = abi^2;
      bj(i) = 4.0 * i * ( i + bt ) * ( i + alf ) * ( i + ab ) ...
        / ( ( abi - 1.0 ) * abi );
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 5 )

    zemu = gamma ( bt + 1.0 );

    for i = 1 : m
      aj(i) = 2.0 * i - 1.0 + bt;
      bj(i) = i * ( i + bt );
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 6 )

    zemu = gamma ( ( bt + 1.0 ) / 2.0 );

    aj(1:m) = 0.0;

    for i = 1 : m
      bj(i) = ( i + bt * mod ( i, 2 ) ) / 2.0;
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 7 )

    ab = bt;
    zemu = 2.0 / ( ab + 1.0 );

    aj(1:m) = 0.0;

    for i = 1 : m
      abi = i + ab * mod(i,2);
      abj = 2 * i + ab;
      bj(i) = abi * abi / ( abj * abj - 1.0 );
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 8 )

    ab = bt + alf;
    zemu = gamma ( bt + 1.0 ) * gamma ( - ( ab + 1.0 ) ) ...
      / gamma ( - alf );
    apone = bt + 1.0;
    aba = ab * apone;
    aj(1) = - apone / ( ab + 2.0 );
    bj(1) = - aj(1) * ( alf + 1.0 ) / ( ab + 2.0 ) / ( ab + 3.0 );
    for i = 2 : m
      abti = ab + 2.0 * i;
      aj(i) = aba + 2.0 * ( ab + i ) * ( i - 1 );
      aj(i) = - aj(i) / abti / ( abti - 2.0 );
    end

    for i = 2 : m - 1
      abti = ab + 2.0 * i;
      bj(i) = i * ( bt + i ) / ( abti - 1.0 ) * ( alf + i ) ...
        / ( abti^2 ) * ( ab + i ) / ( abti + 1.0 );
    end

    bj(m) = 0.0;
    bj(1:m) =  sqrt ( bj(1:m) );

  end

  return
end
function [ d, z ] = imtqlx ( n, d, e, z )

%*****************************************************************************80
%
%% imtqlx() diagonalizes a symmetric tridiagonal matrix.
%
%  Discussion:
%
%    This routine is a slightly modified version of the EISPACK routine to
%    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
%
%    The authors thank the authors of EISPACK for permission to use this
%    routine.
%
%    It has been modified to produce the product Q' * Z, where Z is an input
%    vector and Q is the orthogonal matrix diagonalizing the input matrix.
%    The changes consist (essentialy) of applying the orthogonal transformations
%    directly to Z as they are generated.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 January 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%    Roger Martin, James Wilkinson,
%    The Implicit QL Algorithm,
%    Numerische Mathematik,
%    Volume 12, Number 5, December 1968, pages 377-383.
%
%  Input:
%
%    integer N, the order of the matrix.
%
%    real D(N), the diagonal entries of the matrix.
%
%    real E(N), the subdiagonal entries of the
%    matrix, in entries E(1) through E(N-1). 
%
%    real Z(N), a vector to be operated on.
%
%  Output:
%
%    real D(N), the diagonal entries of the diagonalized matrix.
%
%    real Z(N), the value of Q' * Z, where Q is the matrix that 
%    diagonalizes the input symmetric tridiagonal matrix.
%
  itn = 30;

  prec = eps;

  if ( n == 1 )
    return
  end

  e(n) = 0.0;

  for l = 1 : n

    j = 0;

    while ( true )

      for m = l : n

        if ( m == n )
          break
        end

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) )
          break
        end

      end

      p = d(l);

      if ( m == l )
        break
      end

      if ( j == itn )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'IMTQLX - Fatal error!\n' );
        fprintf ( 1, '  Iteration limit exceeded.\n' );
        error ( 'IMTQLX - Fatal error!' );
      end

      j = j + 1;
      g = ( d(l+1) - p ) / ( 2.0 * e(l) );
      r =  sqrt ( g * g + 1.0 );
      g = d(m) - p + e(l) / ( g + r8_sign ( g ) * abs ( r ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ii = 1 : mml

        i = m - ii;
        f = s * e(i);
        b = c * e(i);

        if ( abs ( f ) >= abs ( g ) )
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e(i+1) = f * r;
          s = 1.0 / r;
          c = c * s;
        else
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e(i+1) = g * r;
          c = 1.0 / r;
          s = s * c;
        end

        g = d(i+1) - p;
        r = ( d(i) - g ) * s + 2.0 * c * b;
        p = s * r;
        d(i+1) = g + p;
        g = c * r - b;
        f = z(i+1);
        z(i+1) = s * z(i) + c * f;
        z(i) = c * z(i) - s * f;

      end

      d(l) = d(l) - p;
      e(l) = g;
      e(m) = 0.0;

    end

  end

  for ii = 2 : n

     i = ii - 1;
     k = i;
     p = d(i);

     for j = ii : n
       if ( d(j) < p )
         k = j;
         p = d(j);
       end
     end

     if ( k ~= i )
       d(k) = d(i);
       d(i) = p;
       p = z(i);
       z(i) = z(k);
       z(k) = p;
     end

  end

  return
end
function parchk ( kind, m, bt, alf )

%*****************************************************************************80
%
%% parchk() checks parameters BETA and ALPHA for classical weight functions.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 January 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Input:
%
%    integer KIND, the rule.
%    1, Legendre,             (a,b)       1.0
%    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
%    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^bt
%    4, Jacobi,               (a,b)       (b-x)^bt*(x-a)^alf
%    5, Generalized Laguerre, (a,inf)     (x-a)^bt*exp(-b*(x-a))
%    6, Generalized Hermite,  (-inf,inf)  |x-a|^bt*exp(-b*(x-a)^2)
%    7, Exponential,          (a,b)       |x-(a+b)/2.0|^bt
%    8, Rational,             (a,inf)     (x-a)^bt*(x+b)^alf
%
%    integer M, the order of the highest moment to
%    be calculated.  This value is only needed when KIND = 8.
%
%    real BETA, ALPHA, the parameters, if required
%    by the value of KIND.
%
  if ( kind <= 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'PARCHK - Fatal error!\n' );
    fprintf ( 1, '  KIND <= 0.\n' );
    error ( 'PARCHK - Fatal error!' );
  end
%
%  Check BETA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
%
  if ( 3 <= kind && bt <= -1.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'PARCHK - Fatal error!\n' );
    fprintf ( 1, '  3 <= KIND and BETA <= -1.\n' );
    error ( 'PARCHK - Fatal error!' );
  end
%
%  Check ALPHA for Jacobi.
%
  if ( kind == 4 && alf <= -1.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'PARCHK - Fatal error!\n' );
    fprintf ( 1, '  KIND == 4 and ALPHA <= -1.0.\n' );
    error ( 'PARCHK - Fatal error!' );
  end
%
%  Check BETA and ALPHA for rational.
%
  if ( kind == 8 )
    tmp = bt + alf + m + 1.0;
    if ( 0.0 <= tmp || tmp <= alf )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'PARCHK - Fatal error!\n' );
      fprintf ( 1, '  KIND == 8 but condition on BETA and ALPHA fails.\n' );
      error ( 'PARCHK - Fatal error!' );
    end
  end

  return
end
function value = r8_sign ( x )

%*****************************************************************************80
%
%% r8_sign() returns the sign of an R8.
%
%  Discussion:
%
%    The value is +1 if the number is positive or zero, and it is -1 otherwise.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    21 March 2004
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real X, the number whose sign is desired.
%
%  Output:
%
%    real VALUE, the sign of X.
%
  if ( 0 <= x )
    value = +1.0;
  else
    value = -1.0;
  end

  return
end
function r8mat_write ( output_filename, m, n, table )

%*****************************************************************************80
%
%% r8mat_write() writes an R8MAT file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 August 2009
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    string OUTPUT_FILENAME, the output filename.
%
%    integer M, the spatial dimension.
%
%    integer N, the number of points.
%
%    real TABLE(M,N), the points.
%

%
%  Open the file.
%
  output_unit = fopen ( output_filename, 'wt' );

  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8MAT_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'R8MAT_WRITE - Error!' );
  end
%
%  Write the data.
%
%  For smaller data files, and less precision, try:
%
%     fprintf ( output_unit, '  %14.6f', table(i,j) );
%
  for j = 1 : n
    for i = 1 : m
      fprintf ( output_unit, '  %24.16f', table(i,j) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  Close the file.
%
  fclose ( output_unit );

  return
end
function rule_write ( order, filename, x, w, r )

%*****************************************************************************80
%
%% rule_write() writes a quadrature rule to a file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    18 February 2010
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer ORDER, the order of the rule.
%
%    string FILENAME, specifies the output files.
%    write files 'filename_w.txt', 'filename_x.txt', 'filename_r.txt' defining 
%    weights, abscissas, and region.
%
%    real X(ORDER), the abscissas.
%
%    real W(ORDER), the weights.
%
%    real R(2), the region.
%
  filename_x = strcat ( filename, '_x.txt' );
  filename_w = strcat ( filename, '_w.txt' );
  filename_r = strcat ( filename, '_r.txt' );

  fprintf ( 1, '\n' );
  fprintf ( 1,'  Creating quadrature files.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  "Root" file name is   "%s".\n', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Weight file will be   "%s".\n', filename_w );
  fprintf ( 1, '  Abscissa file will be "%s".\n', filename_x );
  fprintf ( 1, '  Region file will be   "%s".\n', filename_r );

  r8mat_write ( filename_w, 1, order, w' );
  r8mat_write ( filename_x, 1, order, x' );
  r8mat_write ( filename_r, 1, 2,     r' );

  return
end
function [ t, wts ] = scqf ( nt, t, mlt, wts, nwts, ndx, kind, bt, ...
  alf, a, b )

%*****************************************************************************80
%
%% scqf() scales a quadrature formula to a nonstandard interval.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 February 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Input:
%
%    integer NT, the number of knots.
%
%    real T(NT), the original knots.
%
%    integer MLT(NT), the multiplicity of the knots.
%
%    real WTS(NWTS), the weights.
%
%    integer NWTS, the number of weights.
%
%    integer NDX(NT), used to index the array WTS.
%    For more details see the comments in CAWIQ.
%
%    integer KIND, the rule.
%    1, Legendre,             (a,b)       1.0
%    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
%    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^bt
%    4, Jacobi,               (a,b)       (b-x)^bt*(x-a)^alf
%    5, Generalized Laguerre, (a,+oo)     (x-a)^bt*exp(-b*(x-a))
%    6, Generalized Hermite,  (-oo,+oo)   |x-a|^bt*exp(-b*(x-a)^2)
%    7, Exponential,          (a,b)       |x-(a+b)/2.0|^bt
%    8, Rational,             (a,+oo)     (x-a)^bt*(x+b)^alf
%    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
%
%    real BETA, the value of Alpha, if needed.
%
%    real ALPHA, the value of Beta, if needed.
%
%    real A, B, the interval endpoints.
%
%  Output:
%
%    real T(NT), the scaled knots.
%
%    real WTS(NWTS), the scaled weights.
%
  temp = eps;

  parchk ( kind, 1, bt, alf )

  if ( kind == 1 )

    al = 0.0;
    be = 0.0;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  elseif ( kind == 2 )

    al = -0.5;
    be = -0.5;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  elseif ( kind == 3 )

    al = bt;
    be = bt;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  elseif ( kind == 4 )

    al = bt;
    be = alf;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  elseif ( kind == 5 )

    if ( b <= 0.0 )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  B <= 0.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = a;
    slp = 1.0 / b;
    al = bt;
    be = 0.0;

  elseif ( kind == 6 )

    if ( b <= 0.0 )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  B <= 0.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = a;
    slp = 1.0 / sqrt ( b );
    al = bt;
    be = 0.0;

  elseif ( kind == 7 )

    al = bt;
    be = 0.0;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  elseif ( kind == 8 )

    if ( a + b <= 0.0 )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  A + B <= 0.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = a;
    slp = a + b;
    al = bt;
    be = alf;

  elseif ( kind == 9 )

    al = 0.5;
    be = 0.5;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  end

  p = slp^( al + be + 1.0 );

  for k = 1 : nt

    t(k) = shft + slp * t(k);
    l = abs ( ndx(k) );

    if ( l ~= 0 )
      tmp = p;
      for i = l : l + mlt(k) - 1
        wts(i) = wts(i) * tmp;
        tmp = tmp * slp;
      end
    end

  end

  return
end
function [ t, wts ] = sgqf ( nt, aj, bj, zemu )

%*****************************************************************************80
%
%% sgqf() computes knots and weights of a Gauss Quadrature formula.
%
%  Discussion:
%
%    This routine computes all the knots and weights of a Gauss quadrature
%    formula with simple knots from the Jacobi matrix and the zero-th
%    moment of the weight function, using the Golub-Welsch technique.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 February 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Input:
%
%    integer NT, the number of knots.
%
%    real AJ(NT), the diagonal of the Jacobi matrix.
%
%    real BJ(NT), the subdiagonal of the Jacobi
%    matrix, in entries 1 through NT-1.  On BJ has been overwritten.
%
%    real ZEMU, the zero-th moment of the weight function.
%
%  Output:
%
%    real T(NT), the knots.
%
%    real WTS(NT), the weights.
%

%
%  Exit if the zero-th moment is not positive.
%
  if ( zemu <= 0.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'SGQF - Fatal error!\n' );
    fprintf ( 1, '  ZEMU <= 0.\n' );
    error ( 'SGQF - Fatal error!' );
  end
%
%  Set up vectors for IMTQLX.
%
  wts = zeros ( nt, 1 );
  wts(1) = sqrt ( zemu );
  wts(2:nt) = 0.0;
%
%  Diagonalize the Jacobi matrix.
%
  [ t, wts ] = imtqlx ( nt, aj, bj, wts );

  wts(1:nt) = wts(1:nt).^2;

  return
end
function timestamp ( )

%*****************************************************************************80
%
%% timestamp() prints the current YMDHMS date as a timestamp.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 February 2003
%
%  Author:
%
%    John Burkardt
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end

