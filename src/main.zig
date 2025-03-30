//! By convention, main.zig is where your main function lives in the case that
//! you are building an executable. If you are making a library, the convention
//! is to delete this file and start with root.zig instead.
const std = @import("std");

/// This imports the separate module containing `root.zig`. Take a look in `build.zig` for details.
const lib = @import("kalei-zig_lib");
const Presentation = lib.group.Presentation;
const all_group_elements = lib.group.all_group_elements;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();

    var pres = Presentation.init(allocator, 3);
    defer pres.deinit();

    try pres.add_relation(&[_]usize{0} ** 2);
    try pres.add_relation(&[_]usize{1} ** 2);
    try pres.add_relation(&[_]usize{2} ** 2);
    try pres.add_relation(&[_]usize{ 0, 2 } ** 2);
    try pres.add_relation(&[_]usize{ 0, 1 } ** 4);
    try pres.add_relation(&[_]usize{ 1, 2 } ** 3);

    const elements = try all_group_elements(allocator, pres);
    defer elements.deinit();

    // Prints to stderr (it's a shortcut based on `std.io.getStdErr()`)
    std.debug.print("All your {s} are belong to us.\n", .{"codebase"});

    // stdout is for the actual output of your application, for example if you
    // are implementing gzip, then only the compressed bytes should be sent to
    // stdout, not any debugging messages.
    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("Run `zig build test` to run the tests.\n", .{});

    try bw.flush(); // Don't forget to flush!
}
