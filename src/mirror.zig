const std = @import("std");
const Allocator = std.mem.Allocator;
const ArrayList = std.ArrayList;
const panic = std.debug.panic;
const pi = std.math.pi;
const expect = std.testing.expect;

const matrix = @import("matrix.zig");
const Mat = matrix.Mat;

const CoxeterGroup = struct {
    const Self = @This();

    orders: ArrayList(u8),

    pub fn deinit(self: Self) void {
        self.orders.deinit();
    }

    pub fn order(self: Self, i: usize, j: usize) u8 {
        if (i == j) {
            return 2;
        } else if (i < j) {
            return self.orderUnchecked(i, j);
        } else {
            return self.orderUnchecked(j, i);
        }
    }

    fn orderUnchecked(self: Self, i: usize, j: usize) u8 {
        if (i >= j) panic("Precondition unsatisfied", .{});

        if (i + 1 == j) {
            return self.orders.items[i];
        } else {
            return 2;
        }
    }

    pub fn dihedralAngle(self: Self, i: usize, j: usize) f64 {
        if (i == j) return 0.0;
        const ord: f64 = @floatFromInt(self.order(i, j));
        const angle: f64 = pi / ord;
        return angle;
    }

    pub fn size(self: Self) usize {
        return self.orders.items.len + 1;
    }
};

const Normal = ArrayList(f64);

fn gen_mirrors(allocator: Allocator, group: CoxeterGroup) !ArrayList(Normal) {
    var normals = ArrayList(Normal).init(allocator);
    errdefer {
        for (normals.items) |norm| {
            norm.deinit();
        }
        normals.deinit();
    }

    // This is the very first normal
    try normals.append(Normal.init(allocator));
    try normals.items[0].append(1.0);

    // The other normals are created in order
    for (1..group.size()) |i| {
        // Start by "extending" each existing vector into a new dimension
        for (normals.items) |*norm| {
            try norm.append(0.0);
        }

        // Create the new normal equal to (1, 0, 0, ...) aka the first normal
        const newNorm = try normals.items[0].clone();
        errdefer newNorm.deinit();

        // Rotate this normal to have the correct dihedral angle with each existing normal
        // TODO: Note that this isn't correct for now
        for (0.., normals.items) |j, norm| {
            _ = norm;
            const order: f64 = @floatFromInt(group.order(i, j));
            const dihedralAngle: f64 = pi / order;

            rotate(newNorm, dihedralAngle, 0, 1);
        }

        try normals.append(newNorm);
    }

    return normals;
}

fn rotate(normal: Normal, angle: f64, dim1: usize, dim2: usize) void {
    // Rotation matrix:
    // cos(a), -sin(a)
    // sin(a), cos(a)
    const c = @cos(angle);
    const s = @sin(angle);

    std.debug.print("{} {}, {} {}\n", .{ dim1, dim2, c, s });

    // Only dimensions dim1 and dim2 are affected by the rotation matrix
    const new1 = c * normal.items[dim1] - s * normal.items[dim2];
    normal.items[dim2] = s * normal.items[dim1] + c * normal.items[dim2];
    normal.items[dim1] = new1;
}

fn dotProduct(n1: Normal, n2: Normal) f64 {
    var product: f64 = 0.0;
    for (n1.items, n2.items) |a1, a2| {
        product += a1 * a2;
    }
    return product;
}

fn normalLength(n: Normal) f64 {
    return @sqrt(dotProduct(n, n));
}

fn angleBetween(n1: Normal, n2: Normal) f64 {
    const dot = dotProduct(n1, n2);
    const cos = dot / normalLength(n1) / normalLength(n2);
    return std.math.acos(cos);
}

fn checkAngles(normals: ArrayList(Normal), group: CoxeterGroup) !void {
    try expect(normals.items.len == group.size());

    for (0.., normals.items) |i, norm1| {
        for (0.., normals.items) |j, norm2| {
            const expected = group.dihedralAngle(i, j);
            const actual = angleBetween(norm1, norm2);
            std.debug.print("{any} {any}, {} {}\n", .{ norm1.items, norm2.items, expected, actual });
            try std.testing.expectApproxEqAbs(expected, actual, 0.001);
        }
    }
}

fn printNormals(normals: ArrayList(Normal)) void {
    for (normals.items) |norm| {
        std.debug.print("{any}, {}\n", .{ norm.items, normalLength(norm) });
    }
}

fn addZeroColumn(allocator: Allocator, normals: *Mat(f64)) !void {
    const normalsNew = try Mat(f64).init(allocator, normals.rows, normals.cols + 1);
    errdefer normalsNew.deinit();

    for (0..normals.rows) |i| {
        for (0..normals.cols) |j| {
            normalsNew.getPtr(i, j).* = normals.get(i, j);
        }
    }
    for (0..normals.rows) |i| {
        normalsNew.getPtr(i, normals.rows).* = 0.0;
    }

    normals.*.deinit();

    normals.* = normalsNew;
}

fn appendRow(allocator: Allocator, mat: *Mat(f64), row: Mat(f64)) !void {
    if (row.rows != 1) {
        std.debug.print("{d} {d}\n", .{ row.rows, row.cols });
        return error.NotVector;
    }
    if (mat.cols != row.cols) return error.DimensionMismatch;

    const newMat = try Mat(f64).init(allocator, mat.rows + 1, mat.cols);
    errdefer newMat.deinit();

    for (0..mat.rows) |i| {
        for (0..mat.cols) |j| {
            newMat.getPtr(i, j).* = mat.get(i, j);
        }
    }
    for (0..row.cols) |i| {
        newMat.getPtr(mat.rows, i).* = row.get(0, i);
    }

    mat.deinit();
    mat.* = newMat;
}

fn findNextNormal(allocator: Allocator, pinv: Mat(f64), b: Mat(f64), normals: Mat(f64)) !Mat(f64) {
    const dim = normals.cols;

    const base = try matrix.mul(f64, pinv, b);
    defer base.deinit();

    const tmp1 = try Mat(f64).identity(allocator, dim);
    defer tmp1.deinit();

    const tmp2 = try matrix.mul(f64, pinv, normals);
    defer tmp2.deinit();

    const A = try matrix.sub(f64, tmp1, tmp2);
    defer A.deinit();

    const w = try Mat(f64).init(allocator, dim, 1);
    defer w.deinit();
    w.fill(1.0);

    var scale: f64 = 1.0;

    const Aw = try matrix.mul(f64, A, w);
    defer Aw.deinit();

    for (0..10000) |_| {
        const tmp3 = try matrix.scalarMul(f64, scale, Aw);
        defer tmp3.deinit();

        const x = try matrix.add(f64, base, tmp3);
        defer x.deinit();

        const norm: f64 = try matrix.vectorNorm(f64, x);
        const tmp4 = 4 * (norm - 1);

        var dp: f64 = 0.0;
        for (0..dim) |i| {
            dp += x.get(i, 0) * Aw.get(i, 0);
        }

        const derivative = tmp4 * dp;

        scale -= 0.05 * derivative;
    }

    const tmp3 = try matrix.scalarMul(f64, scale, Aw);
    defer tmp3.deinit();

    const x = try matrix.add(f64, base, tmp3);
    defer x.deinit();

    const norm: f64 = try matrix.vectorNorm(f64, x);

    std.debug.print("{}\n", .{norm});

    return try matrix.scalarMul(f64, 1 / norm, x);

    //     var scaleLower: f64 = 0.0;
    //     var scaleUpper: f64 = 1.0;

    //     // There's an assertion here in the python code
    //     // > # we might need to begin with a higher scale
    //     // > assert np.linalg.norm(base + A @ (scale_upper * w)) > 1.0
    //     // I am ignoring it because my matrix code is annoying to work with

    //     const epsilon = 1e-15;

    //     for (0..100) |_| {
    //         const scale = (scaleLower + scaleUpper) / 2.0;

    //         const tmp3 = try matrix.scalarMul(f64, scale, w);
    //         defer tmp3.deinit();

    //         const tmp4 = try matrix.mul(f64, A, tmp3);
    //         defer tmp4.deinit();

    //         const x = try matrix.add(f64, base, tmp4);
    //         defer x.deinit();

    //         const norm: f64 = try matrix.vectorNorm(f64, x);

    //         std.debug.print("{}\n", .{norm});

    //         if (norm < 1.0 - epsilon) {
    //             scaleLower = scale;
    //         } else if (norm > 1.0 + epsilon) {
    //             scaleUpper = scale;
    //         } else {
    //             return try matrix.scalarMul(f64, 1 / norm, x);
    //         }
    //     }

    //     panic("damn", .{});
}

fn findNormals(allocator: Allocator, group: CoxeterGroup) !void {
    var normals = try Mat(f64).init(allocator, 1, 1);
    defer normals.deinit();

    normals.getPtr(0, 0).* = 1.0;

    while (normals.rows < group.size()) {
        // Add a column of zeros to the normals
        try addZeroColumn(allocator, &normals);

        const numNormals = normals.rows;

        // Set up matrix equation N @ x = b, where:
        // N = normals
        // x = unknown normal
        // b = cosines of dihedral angles

        const b = try Mat(f64).init(allocator, numNormals, 1);
        defer b.deinit();
        for (0..numNormals) |i| {
            const angle = group.dihedralAngle(numNormals - 1, i);
            b.getPtr(i, 0).* = @cos(angle);
        }

        const pinv = try matrix.pinv(normals);
        defer pinv.deinit();

        const x = try findNextNormal(allocator, pinv, b, normals);
        defer x.deinit();

        const xt = try matrix.transpose(f64, x);
        defer xt.deinit();

        try appendRow(allocator, &normals, xt);
    }
}

test {
    const allocator = std.testing.allocator;

    var items = ArrayList(u8).init(allocator);
    try items.append(4);
    try items.append(3);
    const group = CoxeterGroup{
        .orders = items,
    };
    defer group.deinit();

    try findNormals(allocator, group);
    // const normals = try gen_mirrors(allocator, group);
    // defer {
    //     for (normals.items) |norm| {
    //         norm.deinit();
    //     }
    //     normals.deinit();
    // }

    // printNormals(normals);
    // try checkAngles(normals, group);
}
