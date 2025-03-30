//! By convention, main.zig is where your main function lives in the case that
//! you are building an executable. If you are making a library, the convention
//! is to delete this file and start with root.zig instead.
const std = @import("std");

/// This imports the separate module containing `root.zig`. Take a look in `build.zig` for details.
const lib = @import("kalei-zig_lib");
const Presentation = lib.group.Presentation;
const all_group_elements = lib.group.all_group_elements;

pub const std_options = std.Options{
    // Set the log level to info
    .log_level = .info,
};

pub fn main() !void {
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
