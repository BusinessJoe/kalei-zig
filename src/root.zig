//! By convention, root.zig is the root source file when making a library. If
//! you are making an executable, the convention is to delete this file and
//! start with main.zig instead.
const std = @import("std");
const testing = std.testing;
pub const group = @import("group/root.zig");
pub const mirror = @import("mirror.zig");
pub const matrix = @import("matrix.zig");

test {
    _ = group;
    _ = mirror;
    _ = matrix;
}
