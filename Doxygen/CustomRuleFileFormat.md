# Custom Rule File Format

The custom rule functionality allows the creation of a sparse grid using a rule other than the ones implemented in the code. The custom rule is defined via a file with tables that list the levels, number of points per level, exactness of the quadrature at each level, points and their associated weights. Currently, the custom rules work only with global grids and hence the interpolant associated with the rule is a global interpolant using Lagrange polynomials.

The custom rule is defined via custom rule file, with the following format:

* **line 1:** should begin with the string **description:** and it should be followed by a string with a short description of the rule. This string is used only for human readability purposes.
* **line 2:** should begin with the string **levels:** followed by an integer indicating the total number of rule levels defined in the file.

After the description and total number of levels have been defined, the file should contain a sequence of integers describing the number of points and exactness, followed by a sequence of floating point numbers listing the points and weights.

* **integers:** is a sequence of integer pairs where the first integer indicates the number of points for the current level and the second integer indicates the exactness of the rule. For example, the first 3 levels of the Gauss-Legendre rule will be described via the sequence (1, 1 , 2 , 3 , 3 , 5), while the first 3 levels of the Clenshaw-Curtis rule will be described via (1 , 1 , 3 , 3 , 5 , 5).
* **floats:** is a sequence of floating point pairs describing the weights and points. The first number of the pair is the quadrature weight, while the second number if the abscissa. The points associated with the first level are listed in the first pairs. The second set of pairs lists the points associated with the second level and so on.


Here is an example of Gauss-Legendre 3 level rule for reference purposes:
```
description: Gauss-Legendre rule
levels: 3
1 1 2 3 3 5
2.0 0.0
1.0 -0.5774 1.0 0.5774
0.5556 -0.7746 0.8889 0.0 0.5556 0.7746
```

Similarly, a level 3 Clenshaw-Curtis rule can be defined as
```
description: Clenshaw-Curtis rule
levels: 3
1 1 3 3 5 5
2.0 0.0
0.333 1.0 1.333 0.0 0.333 -1.0
0.8 0.0 0.067 -1.0 0.067 1.0 0.533 -0.707 0.533 0.707
```

Several notes on the custom rule file format:
* Tasmanian works with double precision and hence a custom rule should be defined with the corresponding number of significant digits. The examples above are for illustrative purposes only.
* The order of points within each level is irrelevant. Tasmanian will internally index the points.
* Points that are within distance of *1.E-12* of each other will be treated as the same point. Thus, repeated (nested) points can be automatically handled by the code. The tolerance can be adjusted in **tsgMathUtils.hpp** by modifying the **num_tol** constant,
* Naturally, Tasmanian cannot create a sparse grid that requires a one dimensional rule with level higher than what is provided in the file. Predicting the required number of levels can be hard in the case of anisotropic grids, the code will raise a `std::runtime_error` if the custom rule does not provide a sufficient number of points.
* The exactness constants are used only if **qptotal** or **qpcurved** types are used and the indexes of the polynomial space, i.e., **getPolynomialIndexes()**. If quadrature rules are not used, then the exactness integers can be set to 0.
* The quadrature weights are used only if integration is performed. If no quadrature or integration is used, then the weights can all be set to 0.
* If a custom rule is used together with **setDomainTransform()**, then the transform will assume that the rule is defined on the canonical interval `[-1,1]`. A custom rule can be defined on any arbitrary interval, however, for any interval different from `[-1,1]` the **setDomainTransform()** functions should not be used.
* Tasmanian comes with an example custom rule file that defines 9 levels of the Gauss-Legendre-Patterson rule, a.k.a., nested Gauss-Legendre rule.
