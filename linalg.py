class Matrix:
	def __init__(self, rows, cols, entries=None):
		self.rows = rows
		self.cols = cols
		if entries is None:
			self.entries = [0] * (rows * cols)
		elif len(entries) == rows * cols:
			self.entries = entries

	@staticmethod
	def from_formula(rows, cols, f):
		result = Matrix(rows, cols)
		for i in range(rows):
			for j in range(cols):
				result[i, j] = f(i, j)
		return result

	def copy(self):
		return Matrix.from_formula(self.rows, self.cols, lambda i, j: self[i, j])

	def get_entry(self, row, col):
		if 0 <= row < self.rows and 0 <= col < self.cols:
			return self.entries[row * self.cols + col]
		print(self)
		print(row, col)
		raise IndexError

	def __getitem__(self, index):
		return self.get_entry(index[0], index[1])

	def set_entry(self, row, col, value):
		if 0 <= row < self.rows and 0 <= col < self.cols:
			self.entries[row * self.cols + col] = value
		else:
			print(self)
			print(row, col)
			print(self.rows, self.cols)
			raise IndexError

	def __setitem__(self, index, value):
		return self.set_entry(index[0], index[1], value)
	

	def string(self):
		result = []
		for i in range(self.rows):
			for j in range(self.cols):
				result.append(str(self.get_entry(i, j)))
				result.append('\t')
			result.append('\n')
		return "".join(result)

	def __str__(self):
		return self.string()

	def __repr__(self):
		return self.string()

	# LINEAR SYSTEMS AND GAUSSIAN ELIMINATION

	def augment(self, other):
		def f(i, j):
			if j < self.cols:
				return self[i, j]
			return other[i, j - self.cols]
		return Matrix.from_formula(self.rows, self.cols + other.cols, f)

	def slice(self, r0, r1, c0, c1):
		return Matrix.from_formula(r1 - r0, c1 - c0, lambda i, j: self[r0 + i, c0 + j])

	def __or__(self, other):
		return self.augment(other)

	def ero_add(self, from_row, to_row, mult):
		for col in range(self.cols):
			self[to_row, col] += self[from_row, col] * mult

	def ero_mult(self, row, mult):
		for col in range(self.cols):
			self[row, col] *= mult

	def ero_swap(self, row1, row2):
		for col in range(self.cols):
			self[row1, col], self[row2, col] = self[row2, col], self[row1, col]

	def ref(self):
		r0 = c = 0
		while r0 < self.rows and c < self.cols:
			r1 = r0
			while r1 < self.rows and self[r1, c] == 0:
				r1 += 1
			if r1 == self.rows:
				c += 1
				continue
			if r0 != r1:
				self.ero_add(r1, r0, 1)
			for r in range(r1 + 1, self.rows):
				self.ero_add(r0, r, -self[r, c] / self[r0, c])
			r0 += 1
			c += 1

	def get_ref(self):
		result = self.copy()
		result.ref()
		return result

	def rref(self):
		self.ref()
		r = self.rows - 1
		while r >= 0:
			c = 0
			while c < self.cols and self[r, c] == 0:
				c += 1
			if c == self.cols:
				r -= 1
				continue
			self.ero_mult(r, 1 / self[r, c])
			for r1 in range(r):
				self.ero_add(r, r1, -self[r1, c])
			r -= 1

	def get_rref(self):
		result = self.copy()
		result.rref()
		return result

	def get_pivots(self):
		pivots = set()
		c = 0
		for r in range(self.rows):
			while c < self.cols and self[r, c] == 0:
				c += 1
			if c == self.cols:
				break
			pivots.add(c)
		return pivots

	def get_pivot_indices(self):
		pivots = self.get_pivots()
		result = []
		j = 0
		for i in range(self.cols):
			if i in pivots:
				result.append(j)
				j += 1
			else:
				result.append(-1)
		return result

	def get_nonzero_rows(self):
		for row in range(self.rows - 1, -1, -1):
			for col in range(self.cols):
				if self[row, col]:
					return row + 1
		return 0

	def null(self):
		# returns a matrix whose column space is the null space of self
		mat = self.get_rref()
		pivots = mat.get_pivots()
		result_cols = mat.cols - len(pivots)
		result = Matrix(mat.cols, result_cols)
		r = 0
		for c in range(mat.cols):
			if c in pivots:
				c1 = c - r
				for c0 in range(c + 1, mat.cols):
					if c0 in pivots:
						continue
					result[c, c1] = -mat[r, c0]
					c1 += 1
				r += 1
			else:
				result[c, c - r] = 1
		return result

	# MATRICES

	@staticmethod
	def identity(order):
		return Matrix.from_formula(order, order, lambda i, j: 1 if i == j else 0)

	def is_square(self):
		return self.rows == self.cols

	def equal(self, other):
		if self.rows != other.rows and self.cols != other.cols:
			return False
		return self.entries == other.entries

	def __eq__(self, other):
		return self.equal(other)

	def add(self, other):
		if self.rows != other.rows or self.cols != other.cols:
			raise ArithmeticError
		return Matrix.from_formula(self.rows, self.cols, lambda i, j: self[i, j] + other[i, j])

	def __add__(self, other):
		return self.add(other)

	def subtract(self, other):
		if self.rows != other.rows or self.cols != other.cols:
			raise ArithmeticError
		return Matrix.from_formula(self.rows, self.cols, lambda i, j: self[i, j] - other[i, j])

	def __sub__(self, other):
		return self.subtract(other)

	def scale(self, scalar):
		return Matrix.from_formula(self.rows, self.cols, lambda i, j: scalar * self[i, j])

	def mult(self, other):
		if self.cols != other.rows:
			raise ArithmeticError
		return Matrix.from_formula(self.rows, other.cols, lambda i, j: sum(self[i, k] * other[k, j] for k in range(self.cols)))


	def __mul__(self, other):
		if isinstance(other, Matrix):
			return self.mult(other)
		else:
			return self.scale(other)

	def power(self, n):
		if not self.is_square():
			raise ArithmeticError
		result = Matrix.identity(self.rows)
		if n < 0:
			return self.inv().power(-n)
		for i in range(n):
			result *= self
		return result

	def __pow__(self, n):
		return self.power(n)

	def transpose(self):
		return Matrix.from_formula(self.cols, self.rows, lambda i, j: self[j, i])

	def det(self):
		if not self.is_square():
			raise ArithmeticError
		mat = self.get_ref()
		result = 1
		for i in range(self.rows):
			result *= mat[i, i]
		return result


	def inv(self):
		I = Matrix.identity(self.rows)
		aug = self | I
		aug.rref()
		left = aug.slice(0, self.rows, 0, self.rows)
		right = aug.slice(0, self.rows, self.rows, 2 * self.rows)
		if left == I:
			return right
		raise ArithmeticError

	

def solve_matrix_equation(A, B):
	m, p, n = A.rows, A.cols, B.cols

	aug = A | B
	aug.rref()

	# if (A0 | B0) has a pivot entry on the RHS, the matrix equation is inconsistent
	pivots = aug.get_pivots()
	for col in pivots:
		if col >= p:
			return ("inconsistent", None)

	# check for non-pivot column in the LHS
	has_non_pivot = False
	for i in range(p):
		if i not in pivots:
			has_non_pivot = True

	# create the matrices A0, B0
	A0 = aug.slice(0, m, 0, p)
	B0 = aug.slice(0, m, p, p + n)

	# if every column on the LHS is a pivot column, then the matrix equation has a unique solution
	if not has_non_pivot:
		return ("Unique Solution", B0.slice(0, p, 0, n))

	# if the LHS has a non pivot column, then the matrix has infinitely many solutions
	pivot_indices = A0.get_pivot_indices()
	X0 = Matrix.from_formula(p, n, lambda i, j: 0 if pivot_indices[i] == -1 else B0[pivot_indices[i], j])
	V = A0.null()

	return ("Infinitely Many Solutions", X0, V)







A = Matrix(2, 3, [1,-1,1,2,-1,-1])
B = Matrix.identity(2)

print(solve_matrix_equation(A, B))