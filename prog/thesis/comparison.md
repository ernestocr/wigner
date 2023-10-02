# Recollection and comparison of methods.

There doesn't seem to be a unified way of setting up the phase space structure and constructing Wigner functions. Vourdas' methodology is very attractive because it mimics the construction in the continuous case but has an unfortunate disadvantage, he ignores even dimensions from the start, stating technical difficulties:

- In chapter 4 he sets everything up considering odd dimension $d$. He considers the discrete phase space as a cartesian product of the finite ring $\mathbb Z_d$.
- The $x$ axis corresponds to position, which he defines as any orthonormal basis, and calls them *position states* which are labeled by the finite ring elements.
- Then he defines a discrete Fourier transform $F$ using the additive characters $\omega$. Using $F$ he obtains the 'conjugate' variable $p$ from the position states, which he calls the *momentum states*, and are also labelled by the finite ring elements.
- The next step is to define displacement operators. He does this for every point on the phase space and makes no mention of the grouping by translations. The displacement operators are defined using both conjugate variable states, and as in the continuous case, form a special group, which is a discrete version of the Heisenberg-Weyl group.
- After this, he defines the displaced parity operators, and just like the continuous case, he first defines the parity operator around the origin (which is just the Fourier transform squared) and then translates the parity operator using the displacement operators to obtain the parity operator of any point.
- Finally, given the displaced parity operators (which are the point operators that Wootters mentions), he defines the Wigner function just as in the continuous case, i.e., the trace of the product of density and displaced parity opertors.

As we can see, up till now, there is no mention of mutually unbiased basis, all that is required to construct a Wigner function is that $d$ be odd. The methodology mimics the continuous case very closely, Vourdas even defines the Radon transform, and proves that the Wigner function as defined has the marginal properties with respect to the conjugate bases!

So, where does it go wrong for even integers? Since Vourdas offers no analysis, we check via computations what might go wrong, but can also intuit that there are problems with the definitions that use characters. For $d = 2$, the finite Fourier transform is defined appropriately, as it has the properties we wish. The basic displacement operators $X$ and $Z$ are also correctly defined, and as they should, are equal to the Pauli matrices $X$ and $Z$. The first problem occurs in the definition of the general displacement operators, because Vourdas uses the mutliplicative inverse of 2, and the ring $\mathbb Z_2$ doesn't even have the number 2. For $d > 2$ even, the Fourier transform is ill-defined, it is not unitary anymore, $FF^* \neq I$ and doesn't satisfy $FF^4 = I$. The $Z$ operator is well defined but $X$ is not, it is not the basis shift operator anymore. This makes sense because the $X$ operator is defined using the momentum basis which in turns requires the additive characters. But wait, for $d = 9$ it doesn't work, is it because I'm not using rings but fields? Wootters ignores the phase factors of the general displacement operators, but Vourdas includes them directly:
$$D(\alpha,\beta) = Z^\alpha X^\beta \omega(-2^{-1}\alpha\beta), \quad \alpha,\beta \in \mathbb Z_d.$$
Does it have to do with the commutation relations? He does mention that the operators $D(\alpha,\beta)\omega(\gamma)$ form a representation of the $HW[\mathbb Z_d]$ Heisenberg-Weyl group...

Before introducing the Wigner and Weyl functions he then proceeds to talk about some special groups, which he doesn't really use in the chapter. He defines a unitary representation of the $Sp[2,\mathbb Z_d]$ group, which are the matrices $S(\kappa, \lambda | \mu, \nu)$. And after this he actually combines both the displacement operators and the symplectic transformations, to form a so denoted $HWSp[\mathbb Z_d]$ group of displacements and symplectic transformations. It turns out that this new group is the semi-direct product of the $HW[\mathbb Z_d]$ by the $Sp[2,\mathbb Z_d]$. In relation to MUBs, in chapter 5, Vourdas considers prime dimension systems, which allows us to work with what he calls near-linear finite geometry $\mathbb Z_d \times \mathbb Z_d$. By introducing lines and symplectic transformations of the phase space, he can construct the mutually unbiased bases *directly*, while keeping a correspondence between lines in the phase space and elements of the mutually unbiased bases. The correspondence is simply
$$\mathcal L(\nu) \leftrightarrow \{|\mathbb X(\nu); m\rangle\}, \quad \nu = -1,...,p-1,$$
where $\mathcal L(\nu) = L(1,\nu)$, and $L(\rho,\sigma) = \{(\tau \rho, \tau \sigma) : \tau \in \mathbb Z_d\}$ is a line through the origin.
 
## Composite systems

In chapter 9, Vourdas constructs the discrete phase space system for prime power dimension Hilbert spaces. The setup uses the Galois field as usual and defines the characters using the field trace. Again he does not consider the case where the Galois field is of characteristic 2. Wootters method doesn't care about the characteristic, so we must analyze where it is different. To do this first we mention Vourdas' method for primer power dimension:

- To start off, consider a quantum system with variables in $GF(p^e)$ where $p$ is prime number $p\neq 2$. The elements of $GF(p^e)$ are vectors in $[\mathbb Z_p]^e$ w.r.t. addition, obviously the quantum Galois system has more structure than a simple product of prime quantum systems.
- Vourdas first sets up the Fourier transform like he did before. But instead of the additive characters previously defined he uses new ones, which use the field trace map.
- Momentum states are again obtained by applying the Fourier transform to the position states, which again are a selected orthonormal basis.
- Decomposing any element of the Galois field in terms of a field basis $\{\epsilon^i\}$, $m = \sum m_i \epsilon^i$, Vourdas shows a natural correspondence between the Galois quantum system and product of prime quantum systems:
  $$H[GF(p^e)] \leftrightarrow H[\mathbb Z_p] \otimes \cdots \otimes H[\mathbb Z_p]$$
  where
  $$|X; m\rangle \mapsto |X; m_0\rangle \otimes \cdots \otimes |X; m_{e-1}\rangle.$$

- This decomposition allows us to express the Fourier transform as a sum of position state tensor products, and then define the momentum states as tensor products of momentum states in the prime systems. What is interesting is that in this definition the dual components of a Galois element appear.
- He then defines the basic displacement operators $Z,X$, using the additive character $\chi$ and the Fourier transform.
- Using the tensor expression of the position and momentum states we obtain a tensor expression for the displacement operators.
- Again he defines general displacement operators as $D(\alpha,\beta) = Z(\alpha)X(\beta)\chi(-2^{-1}\alpha\beta)$. Their action on position and momentum states is the same as the prime case.
- Once more, he expresses the displacement operators as tensor products of prime dimension displacement operators, noting that both the basis and dual components appear:
  $$D(\alpha,\beta) = D(\hat{\alpha}_0,\beta_0) \otimes \cdots \otimes D(\hat{\alpha}_{e-1},\beta_{e-1}).$$

- Next he defines the symplectic transformations as $p^e \times p^e$ unitary matrices by their action on the $X(\beta)$ and $Z(\alpha)$ operators. These form a unitary representation of the $Sp[2, GF(p^e)]$ group.
- An analogous expression for the symplectic transformations in terms of the position states is given using Gauss sums (modified to use the trace).
- Next he defines the partiy operators much in the same way as before. We have a tensor product expression given by the prime dimension parity operators, where once again $P(0,0) = F^2$.
- He doesn't even define the Wigner function because it is just the same as before using the displaced parity operators.
- He ends the section by constructing the $p^e + 1$ mutually unbiased bases in the same way as before, using the symplectic transformations $S(0,-1|1,\nu)$. He notes that the geometry $GF(p^e) \times GF(p^e)$ has $p^e+1$ lines through the origin, and that there is a duality between the $p^e+1$ bases and the $p^e+1$ lines through the origin where $$\mathcal B(\nu) \leftrightarrow \mathcal L(\nu).$$

So how does this compare to Wootters' method? In the 2004 paper, Wootters starts directly from a quantum system of dimension $p^e$.

- Wootters starts out with the discrete phase space geometry given by $GF(p^e) \times GF(p^e)$. He notes that there are exactly $N(N+1)$ in the phase space, and each can be grouped into $N+1$ sets of parallel lines, which he calls a striation of the phase space. 
- He then defines geometrical translations of the phase space as $\mathcal T_\alpha \beta = \alpha + \beta$. Of course these translations aren't always intuitive because of the field extension operations.
- He then defines the field basis, field trace and dual basis.
- Wootters method consists of supplying the discrete phase space with a quantum interpretation. He does this by assigning to each line in phase space a specific pure quantum state. This is done by the definition of a so called *quantum net* $Q$. The choice of quantum net is restricted by a few properties. Mainly the property of translational covariance.
- Next he defines translation operators (Vourdas and the rest's displacement operators) corresponding to each phase space translation, he does this by requiring that the composition of translation operators mimics the composition of translations.
- The key to connecting translation vectors $(x,y)$ to individual particles is done by expanding $x$ and $y$ in field bases, using *two* field basis. 
- He chooses the standard basis and first defines the *unit horizontal translation* $(1,0)$ operator as the basis element shifting operator $X$. 
- The *unit vertical translation* $(0,1)$ is then defined as the phase operator.
- Wootters notes that these operators are the *generalized* Pauli matrices.
- Next he defines the rest of unitary translation operators as a tensor product of these basic translation operators:
  $$T(x,y) = X^{x_1} Z^{y_1} \otimes \cdots \otimes X^{x_n} Z^{y_n}.$$
  Notice how Wootters doesn't care about the phase. This is probably what allows him to work with characteristic two.
- Next he restricts the quantum net to satisfy translational covariance, that is to say translations of lines correspond to translations of quantum states.
- This requirement is strong, because for translations that leave a certain line invariant, it amounts to a requirement of commutation between the translation operators.
- The commutativity of the translation operators corresponding to same-invariance translations, along with the fact that they are traceless and trace-orthogonal, implies the definition of a unique basis of simultaneous eigenvecetors up to phase factors.
- Wootters then uses a result from Bandyopadhyay, stating that bases obtained in such a manner are mutually unbiased.
- He then picks out the basis for the basic translation operators. For $T(0,s) = Z^{s_1} \otimes \cdots \otimes Z^{s_n}$, we have the standard basis. For the $T(s,0) = X^{s_1} \otimes \cdots \otimes X^{s_n}$, he obtains the momentum states, (which he doesn't actully mention):
  $$|j\rangle = \frac{1}{\sqrt p}  \sum_{k=1}^p \omega^{jk} |k\rangle.$$

- After obtaining the mutually unbiased bases, which are in correspondence with the striations. He then assigns a specific state $Q(\lambda)$ to each line of each striation. He again asks for translation covariance, which completely determines the quantum net once we *choose* a state to be assigned to each *ray* of each striation.
- It seems Wootters really went backwards on everything. To define the Wigner function, he starts from the property that the Wigner function summed over the line $\lambda$ is the probability that quantum system will be found in the state $Q(\lambda)$.
- This condition completely determines the Wigner function. He uses the same equation as always $W_\alpha = \frac{1}{N} \text{Tr}(\rho A_\alpha)$ where $A_\alpha$ is the point operator (the displaced parity operator in our language), given by 
  $$A_\alpha = \sum_{\lambda \ni \alpha} Q(\lambda) - I.$$
  That is to say, the point operators are defined as the sum of the quantum states of the lines that pass through the point $\alpha$.

## Comparison

The big question is how different are these two methods? Clearly Vourdas' methodology is more natural and mathematically sound, but it seems Wootters is even more general. In Vourdas' method, we don't actually need to define the MUBs to define the Wigner function. In Wootters' method, it is necessary, and there is possibility for different Wigner functions. The fact that Wootters' method doesn't discriminate characteristic two systems is an advantage, but then again having to construct the basis by simultaneous diagonalization is terrible. Godsil and Roy proved that Wootters MUBs and other explicit constructions are unitarily equivalent, I suspect Vourdas' MUBs are also equivalent. My biggset question so far is, how would we go about defining the Wigner function with different MUBs? Wootters' method shows us a way, by creating a quantum net, but, how do we preserve translational covariance when the distinct MUBs may not be eigenbasis of the translation operators? Or must they be if they are unitarily equivalent? What happens if we use not unitarily equivalent MUBs?

To answer some of these questions I'm going to use Vourdas' and Wootters' method and compare step by step.

- First, why can't Vourdas' method work for characteristic two!? 
- Vourdas' tensor expressions for the displacement operators match in the case of prime dimension.
- How does the phase factor $\chi(\alpha\beta)$ affect the displacement operators? We notice that both methods define the displacement operators as a tensor product of powers of basic displacement operators. In Vourdas' method, the $Z$ power correspond to the dual components of $\alpha$ and the $X$ powers correspond to the components of $\beta$. Wooters uses the components of $\alpha$ for the powers of $X$ and a multiple of the dual components for the powers of $Z$. Notice that they both define the order of the product $X$ with $Z$ in opposite fashion. 
- How do the MUBs compare?
- Are they unitarily equivalent?
- How do the displaced parity operators compare to the point operators?
- What happens if we use unitarily equivalent MUBs? Do the marginal properties still hold?

## Klimov's version version for $2^n$ systems.

They note that Wootters version of the discrete Wigner function is done in the same vein as Schwinger, who was the first to realize that the expansion of arbitrary operators in terms of certain operator basis was the crucial concept setting a proper phase-space description. He identified the finite counterpart of the Weyl-Heisenberg group which describes the conjugacy of position and momentum and is used to define a proper $d \times d$ phase-space.