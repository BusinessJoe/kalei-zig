const std = @import("std");
const Allocator = std.mem.Allocator;

pub fn Mat(T: anytype) type {
    return struct {
        const Self = @This();

        rows: usize,
        cols: usize,
        data: []T,
        allocator: Allocator,

        pub fn init(allocator: Allocator, rows: usize, cols: usize) Allocator.Error!Self {
            return Self{
                .rows = rows,
                .cols = cols,
                .data = try allocator.alloc(T, rows * cols),
                .allocator = allocator,
            };
        }

        pub fn deinit(self: Self) void {
            self.allocator.free(self.data);
        }

        pub fn get(self: Self, r: usize, c: usize) T {
            if (r >= self.rows) @panic("row index too large");
            if (c >= self.cols) @panic("col index too large");

            return self.data[r * self.cols + c];
        }

        pub fn getPtr(self: Self, r: usize, c: usize) *T {
            if (r >= self.rows) @panic("row index too large");
            if (c >= self.cols) @panic("col index too large");

            return &self.data[r * self.cols + c];
        }

        pub fn fill(self: Self, value: T) void {
            for (0..self.rows) |i| {
                for (0..self.cols) |j| {
                    self.getPtr(i, j).* = value;
                }
            }
        }

        pub fn clone(self: Self) Allocator.Error!Self {
            return Self{
                .rows = self.rows,
                .cols = self.cols,
                .data = try self.allocator.dupe(T, self.data),
                .allocator = self.allocator,
            };
        }

        pub fn identity(allocator: Allocator, dim: usize) Allocator.Error!Self {
            const mat = try Self.init(allocator, dim, dim);
            errdefer mat.deinit();

            mat.fill(0.0);

            for (0..dim) |i| {
                mat.getPtr(i, i).* = 1.0;
            }

            return mat;
        }
    };
}

fn zero(comptime T: type) T {
    if (T == f32) return 0.0;
    if (T == f64) return 0.0;
    return 0;
}

pub fn add(comptime T: type, a: Mat(T), b: Mat(T)) (error{DimensionMismatch} || Allocator.Error)!Mat(T) {
    if (a.rows != b.rows or a.cols != b.cols) return error.DimensionMismatch;

    var c = try a.clone();
    errdefer c.deinit();

    for (0..c.rows) |i| {
        for (0..c.cols) |j| {
            c.getPtr(i, j).* += b.get(i, j);
        }
    }

    return c;
}

pub fn sub(comptime T: type, a: Mat(T), b: Mat(T)) (error{DimensionMismatch} || Allocator.Error)!Mat(T) {
    if (a.rows != b.rows or a.cols != b.cols) return error.DimensionMismatch;

    var c = try a.clone();
    errdefer c.deinit();

    for (0..c.rows) |i| {
        for (0..c.cols) |j| {
            c.getPtr(i, j).* -= b.get(i, j);
        }
    }

    return c;
}

pub fn scalarMul(comptime T: type, scalar: T, a: Mat(T)) Allocator.Error!Mat(T) {
    var mat = try a.clone();
    errdefer mat.deinit();

    for (0..a.rows) |i| {
        for (0..a.cols) |j| {
            mat.getPtr(i, j).* *= scalar;
        }
    }

    return mat;
}

pub fn mul(comptime T: type, a: Mat(T), b: Mat(T)) (error{DimensionMismatch} || Allocator.Error)!Mat(T) {
    const allocator = a.allocator;

    if (a.cols != b.rows) return error.DimensionMismatch;

    const rows = a.rows;
    const cols = b.cols;

    var c = try Mat(T).init(allocator, rows, cols);
    errdefer c.deinit();

    for (0..rows) |i| {
        for (0..cols) |j| {
            c.getPtr(i, j).* = zero(T);
            for (0..a.cols) |k| {
                c.getPtr(i, j).* += a.get(i, k) * b.get(k, j);
            }
        }
    }

    return c;
}

pub fn transpose(comptime T: type, a: Mat(T)) Allocator.Error!Mat(T) {
    const allocator = a.allocator;

    const rows = a.cols;
    const cols = a.rows;

    var mat = try Mat(T).init(allocator, rows, cols);
    errdefer mat.deinit();

    for (0..rows) |i| {
        for (0..cols) |j| {
            mat.getPtr(i, j).* = a.get(j, i);
        }
    }

    return mat;
}

fn panicDimensionMismatch() noreturn {
    @panic("dimension mismatch!");
}

/// Computes the Moore-Penrose inverse.
pub fn pinv(a: Mat(f64)) Allocator.Error!Mat(f64) {
    const u = 0.0001;

    const at = try transpose(f64, a);
    defer at.deinit();

    var xprev = try at.clone();
    defer xprev.deinit();

    var xnew = try at.clone();
    errdefer xnew.deinit();

    for (0..100) |_| {
        const tmp1 = mul(f64, a, xprev) catch |err| switch (err) {
            error.DimensionMismatch => panicDimensionMismatch(),
            error.OutOfMemory => return error.OutOfMemory,
        };
        defer tmp1.deinit();

        const tmp2 = mul(f64, tmp1, a) catch |err| switch (err) {
            error.DimensionMismatch => panicDimensionMismatch(),
            error.OutOfMemory => return error.OutOfMemory,
        };
        defer tmp2.deinit();

        const tmp3 = sub(f64, a, tmp2) catch |err| switch (err) {
            error.DimensionMismatch => panicDimensionMismatch(),
            error.OutOfMemory => return error.OutOfMemory,
        };
        defer tmp3.deinit();

        const tmp4 = try scalarMul(f64, u, at);
        defer tmp4.deinit();

        const tmp5 = mul(f64, tmp4, tmp3) catch |err| switch (err) {
            error.DimensionMismatch => panicDimensionMismatch(),
            error.OutOfMemory => return error.OutOfMemory,
        };
        defer tmp5.deinit();

        const tmp6 = mul(f64, tmp5, at) catch |err| switch (err) {
            error.DimensionMismatch => panicDimensionMismatch(),
            error.OutOfMemory => return error.OutOfMemory,
        };
        defer tmp6.deinit();

        const tmp7 = add(f64, xprev, tmp6) catch |err| switch (err) {
            error.DimensionMismatch => panicDimensionMismatch(),
            error.OutOfMemory => return error.OutOfMemory,
        };
        defer tmp7.deinit();

        xnew.deinit();
        xnew = try tmp7.clone();

        xprev.deinit();
        xprev = try xnew.clone();
    }

    return xnew;
}

fn rowVectorNorm(comptime T: type, vector: Mat(T)) T {
    var sum = zero(T);
    for (0..vector.cols) |i| {
        const x = vector.get(0, i);
        sum += x * x;
    }
    return @sqrt(sum);
}

fn colVectorNorm(comptime T: type, vector: Mat(T)) T {
    var sum = zero(T);
    for (0..vector.rows) |i| {
        const x = vector.get(i, 0);
        sum += x * x;
    }
    return @sqrt(sum);
}

pub fn vectorNorm(comptime T: type, vector: Mat(T)) error{NotVector}!T {
    if (vector.rows == 1) {
        return rowVectorNorm(T, vector);
    } else if (vector.cols == 1) {
        return colVectorNorm(T, vector);
    } else {
        return error.NotVector;
    }
}

/// Perform the Cholesky decomposition of a matrix, returning the lower-triangular part.
pub fn cholesky_decomp(allocator: Allocator, matrix: Mat(f64)) Allocator.Error!Mat(f64) {
    if (matrix.rows != matrix.cols) {
        @panic("matrix must be square");
    }

    // Implementation of the Cholesky–Banachiewicz algorithm.

    const n = matrix.rows;

    var lower = try Mat(f64).init(allocator, n, n);
    errdefer lower.deinit();

    for (0..n) |i| {
        for (0..i+1) |j| {
            var sum: f64 = 0.0;
            for (0..j) |k| {
                sum += lower.get(i, k) * lower.get(j, k);
            }

            if (i == j) {
                lower.getPtr(i, j).* = @sqrt(matrix.get(i, i) - sum);
            } else {
                lower.getPtr(i, j).* = (1.0 / lower.get(j, j) * (matrix.get(i, j) - sum));
            }
        }
    }

    return lower;
}

test "matrix multiplication int" {
    const allocator = std.testing.allocator;
    const expectEqual = std.testing.expectEqual;

    const a = try Mat(i32).init(allocator, 2, 3);
    defer a.deinit();
    a.fill(0);
    a.getPtr(0, 0).* = 3;
    a.getPtr(0, 1).* = -4;
    a.getPtr(0, 2).* = -9;
    a.getPtr(1, 0).* = -1;
    a.getPtr(1, 1).* = 0;
    a.getPtr(1, 2).* = 8;

    const b = try Mat(i32).init(allocator, 3, 2);
    defer b.deinit();
    b.fill(0);
    b.getPtr(0, 0).* = -2;
    b.getPtr(0, 1).* = 7;
    b.getPtr(1, 0).* = 5;
    b.getPtr(1, 1).* = -6;
    b.getPtr(2, 0).* = 4;
    b.getPtr(2, 1).* = 0;

    const c = try mul(i32, a, b);
    defer c.deinit();

    try expectEqual(2, c.rows);
    try expectEqual(2, c.cols);
    try expectEqual(-62, c.get(0, 0));
    try expectEqual(45, c.get(0, 1));
    try expectEqual(34, c.get(1, 0));
    try expectEqual(-7, c.get(1, 1));
}

test "matrix multiplication float" {
    const allocator = std.testing.allocator;
    const expectEqual = std.testing.expectEqual;
    const expectApproxEqAbs = std.testing.expectApproxEqAbs;

    const a = try Mat(f64).init(allocator, 2, 3);
    defer a.deinit();
    a.fill(0.0);
    a.getPtr(0, 0).* = 3.0;
    a.getPtr(0, 1).* = -4.0;
    a.getPtr(0, 2).* = -9.0;
    a.getPtr(1, 0).* = -1.0;
    a.getPtr(1, 1).* = 0.0;
    a.getPtr(1, 2).* = 8.0;

    const b = try Mat(f64).init(allocator, 3, 2);
    defer b.deinit();
    b.fill(0.0);
    b.getPtr(0, 0).* = -2.0;
    b.getPtr(0, 1).* = 7.0;
    b.getPtr(1, 0).* = 5.0;
    b.getPtr(1, 1).* = -6.0;
    b.getPtr(2, 0).* = 4.0;
    b.getPtr(2, 1).* = 0.0;

    const c = try mul(f64, a, b);
    defer c.deinit();

    try expectEqual(2, c.rows);
    try expectEqual(2, c.cols);
    try expectApproxEqAbs(-62.0, c.get(0, 0), 0.000001);
    try expectApproxEqAbs(45.0, c.get(0, 1), 0.000001);
    try expectApproxEqAbs(34.0, c.get(1, 0), 0.000001);
    try expectApproxEqAbs(-7.0, c.get(1, 1), 0.000001);
}

test "matrix transpose" {
    const allocator = std.testing.allocator;
    const expectEqual = std.testing.expectEqual;

    const a = try Mat(i32).init(allocator, 2, 3);
    defer a.deinit();
    a.fill(0);
    a.getPtr(0, 0).* = 3;
    a.getPtr(0, 1).* = -4;
    a.getPtr(0, 2).* = -9;
    a.getPtr(1, 0).* = -1;
    a.getPtr(1, 1).* = 0;
    a.getPtr(1, 2).* = 8;

    const b = try transpose(i32, a);
    defer b.deinit();
    try expectEqual(3, b.rows);
    try expectEqual(2, b.cols);
    try expectEqual(3, b.get(0, 0));
    try expectEqual(-1, b.get(0, 1));
    try expectEqual(-4, b.get(1, 0));
    try expectEqual(0, b.get(1, 1));
    try expectEqual(-9, b.get(2, 0));
    try expectEqual(8, b.get(2, 1));
}

test "matrix pinv" {
    const allocator = std.testing.allocator;
    const expectEqual = std.testing.expectEqual;
    const expectApproxEqAbs = std.testing.expectApproxEqAbs;

    const a = try Mat(f64).init(allocator, 2, 2);
    defer a.deinit();
    a.fill(0.0);
    a.getPtr(0, 0).* = 3.0;
    a.getPtr(0, 1).* = 2.0;
    a.getPtr(1, 0).* = 6.0;
    a.getPtr(1, 1).* = 4.0;

    const p = try pinv(a);
    defer p.deinit();
    try expectEqual(2, p.rows);
    try expectEqual(2, p.cols);

    const epsilon = 0.000001;
    try expectApproxEqAbs(0.04615385, p.get(0, 0), epsilon);
    try expectApproxEqAbs(0.09230769, p.get(0, 1), epsilon);
    try expectApproxEqAbs(0.03076923, p.get(1, 0), epsilon);
    try expectApproxEqAbs(0.06153846, p.get(1, 1), epsilon);
}
