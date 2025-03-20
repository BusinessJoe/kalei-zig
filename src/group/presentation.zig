const std = @import("std");
const Allocator = std.mem.Allocator;
const ArrayList = std.ArrayList;

/// A finite presentation of a group.
pub const Presentation = struct {
    const Self = @This();

    /// Number of non-identity generators of the group.
    num_gens: usize,
    /// Relations between generators of the group.
    rels: ArrayList(Relation),

    /// Initializes a presentation with the given number of generators and zero relations.
    /// Relations should be added afterwards with `add_relation`.
    /// Deinitialize with `deinit`.
    pub fn init(allocator: Allocator, num_gens: usize) Self {
        return Self{
            .num_gens = num_gens,
            .rels = ArrayList(Relation).init(allocator),
        };
    }

    /// Release all allocated memory.
    pub fn deinit(self: Self) void {
        self.rels.deinit();
    }

    /// Add a relation to the presentation.
    pub fn add_relation(self: *Self, rel: Relation) !void {
        try self.rels.append(rel);
    }
};

/// A relation is a sequence of generators whose product is the identity.
/// Here they are represented by integers which are less than `Presentation.num_gens`.
pub const Relation = []const usize;

// Tests
const test_allocator = std.testing.allocator;

test "use presentation" {
    var pres = Presentation.init(test_allocator, 2);
    defer pres.deinit();
    try pres.add_relation(&[_]usize{ 0, 0 });
    try pres.add_relation(&[_]usize{ 1, 1 });
    try pres.add_relation(&[_]usize{ 0, 1, 0, 1, 0, 1 });
}
