//! By convention, main.zig is where your main function lives in the case that
//! you are building an executable. If you are making a library, the convention
//! is to delete this file and start with root.zig instead.
const std = @import("std");

/// This imports the separate module containing `root.zig`. Take a look in `build.zig` for details.
const lib = @import("kalei-zig_lib");
const Presentation = lib.group.Presentation;
const all_group_elements = lib.group.all_group_elements;

const lapack = @cImport(@cInclude("lapacke.h"));

pub const std_options = std.Options{
    // Set the log level to info
    .log_level = .info,
};

const laint = lapack.lapack_int;

const n: laint = 4;
const nrhs: laint = 1;
const lda: laint = n;
const ldb: laint = n;

const matrix_t = [n * n]f32;
const vector_t = [n]f32;

fn lapack_demo() !void {
    var mat: matrix_t = matrix_t{ //
        1.80, 2.88, 2.05, -0.89, //
        5.25, -2.95, -0.95, -3.80, //
        1.58, -2.69, -2.90, -1.04, //
        -1.11, -0.66, -0.59, 0.80, //
    };  
    var vec: vector_t = vector_t{ 9.52, 24.35, 0.77, -6.22 };
    std.debug.print("{any}\n", .{mat});
    std.debug.print("{any}\n", .{vec});

    var ipiv = [n]laint{ 0, 0, 0, 0 };

    const info = lapack.LAPACKE_sgesv( //
        lapack.LAPACK_COL_MAJOR, // matrix orientation
        n, // number of equations in matrix
        nrhs, // number of colums in result vector
        &mat[0], // the matrix to be diagonalized
        lda, // leading dimension of matrix
        &ipiv[0], // pivot indices
        &vec[0], // result vector
        ldb); // leading dimension of result vector

    std.debug.print("{d}\n", .{info});
    std.debug.print("{any}\n", .{mat});
    std.debug.print("{any}\n", .{ipiv});
    std.debug.print("{any}\n", .{vec});
}

pub fn main() !void {
    try lapack_demo();

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();

    var pres = Presentation.init(allocator, 4);
    defer pres.deinit();

    try pres.add_relation(&[_]usize{0} ** 2);
    try pres.add_relation(&[_]usize{ 0, 1 } ** 5);
    try pres.add_relation(&[_]usize{ 0, 2 } ** 2);
    try pres.add_relation(&[_]usize{ 0, 3 } ** 2);
    try pres.add_relation(&[_]usize{1} ** 2);
    try pres.add_relation(&[_]usize{ 1, 2 } ** 3);
    try pres.add_relation(&[_]usize{ 1, 3 } ** 2);
    try pres.add_relation(&[_]usize{2} ** 2);
    try pres.add_relation(&[_]usize{ 2, 3 } ** 3);
    try pres.add_relation(&[_]usize{3} ** 2);

    const elements = try all_group_elements(allocator, pres);
    defer elements.deinit();

    std.debug.print("{}\n", .{elements.items.len});

    if (elements.items.len != 14400) return error.WrongGroupOrder;
}
