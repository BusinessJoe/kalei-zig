//! Implementation of Todd-Coxeter algorithm to find unique elements of a
//! group given by a finite presentation (`all_group_elements`).

const std = @import("std");
const Allocator = std.mem.Allocator;
const ArrayList = std.ArrayList;
const panic = std.debug.panic;
const log = std.log;
const presentation = @import("presentation.zig");
const Presentation = presentation.Presentation;
const Relation = presentation.Relation;

/// Find all the distinct elements of the group represented
/// by the given presentation. All generators must be involutions.
/// The return elements are represented by a word made from the presentation's generators.
pub fn all_group_elements(allocator: Allocator, pres: Presentation) ![]GeneratorWord {
    var coset_table = CosetTable.init(allocator, pres.num_gens);
    defer coset_table.deinit();

    var rel_tables = ArrayList(RelTable).init(allocator);
    defer {
        for (rel_tables.items) |rel_table| {
            rel_table.deinit();
        }
        rel_tables.deinit();
    }

    for (pres.rels.items) |rel| {
        try rel_tables.append(RelTable.init(allocator, rel));
    }

    // Add first coset to coset_table and rel_tables
    try coset_table.new_coset(1);
    for (rel_tables.items) |*rel_table| {
        try rel_table.new_coset(1);
    }

    // Deductions which have yet to be procesed are stored in a stack
    var deductions = ArrayList(Deduction).init(allocator);
    defer deductions.deinit();

    // TODO: special case if any of the rel_tables have a width (relation length) of 1.
    // The row is already completed and so we have a new deduction.

    var num_cosets: usize = 1;
    while (true) {
        // handle any deductions we have
        while (deductions.pop()) |deduction| {
            std.debug.print("Adding deduction {any}\n", .{deduction});
            if (coset_table.add_deduction(deduction)) |coincidence| {
                // If filling the coset table produces a coincidence Ci = Cj, i < j
                // we need to replace all the instances of j with i.
                // Then fill in all possible entries of the tables

                // Right now we just panic lol
                _ = coincidence;
                panic("Coincidences not implemented", .{});
            }

            for (rel_tables.items) |*rel_table| {
                try rel_table.update_with_coset_table(coset_table, &deductions);
            }
        }

        // the algorithm terminates when all the tables are full
        var all_full = true;
        all_full = all_full and coset_table.is_full();
        for (rel_tables.items) |rel_table| {
            all_full = all_full and rel_table.is_full();
        }
        if (all_full) {
            break;
        }

        // otherwise we add a new coset to the coset table, creating a deduction
        num_cosets += 1;
        std.debug.print("Adding coset {}\n", .{num_cosets});
        try coset_table.new_coset(num_cosets);
        for (rel_tables.items) |*rel_table| {
            try rel_table.new_coset(num_cosets);
        }
        const unknown = coset_table.find_unknown();
        try deductions.append(unknown.to_deduction(num_cosets));
    }

    // Now that the coset table is complete, we can build the group's elements
    return try elements_from_coset_table(allocator, coset_table);
}

fn elements_from_coset_table(allocator: Allocator, coset_table: CosetTable) ![]GeneratorWord {
    var elements = try allocator.alloc(ArrayList(usize), coset_table.num_cosets);
    for (0..coset_table.num_cosets) |i| {
        elements[i] = ArrayList(usize).init(allocator);
    }
    defer allocator.free(elements);
    errdefer {
        for (elements) |word| {
            word.deinit();
        }
    }

    // All elements are initialized to the identity element (an empty list)

    var coset_stack = ArrayList(usize).init(allocator);
    defer coset_stack.deinit();
    try coset_stack.append(1); // We've only added the first coset (the identity) so far

    while (coset_stack.pop()) |coset_in| {
        for (0..coset_table.num_gens) |gen| {
            // The table is filled so this lookup always succeeds
            const coset_out = coset_table.lookup(coset_in, gen).?;
            if (coset_out != 1 and elements[coset_out - 1].items.len == 0) {
                var word = try elements[coset_in - 1].clone();
                try word.append(gen);
                elements[coset_out - 1] = word;
                try coset_stack.append(coset_out);
            }
        }
    }

    var non_opt_elements = try allocator.alloc(GeneratorWord, coset_table.num_cosets);
    for (0..coset_table.num_cosets) |i| {
        non_opt_elements[i] = try elements[i].toOwnedSlice();
    }
    return non_opt_elements;
}

pub const GeneratorWord = []usize;

const CosetTable = struct {
    const Self = @This();

    num_gens: usize,
    num_cosets: usize,
    table_raw: ArrayList(usize),

    pub fn init(allocator: Allocator, num_gens: usize) Self {
        return Self{
            .num_gens = num_gens,
            .num_cosets = 0,
            .table_raw = ArrayList(usize).init(allocator),
        };
    }

    pub fn deinit(self: Self) void {
        self.table_raw.deinit();
    }

    pub fn new_coset(self: *Self, id: usize) !void {
        _ = id;

        for (0..self.num_gens) |_| {
            try self.table_raw.append(0);
        }

        self.num_cosets += 1;
    }

    pub fn add_deduction(self: *Self, deduction: Deduction) ?Coincidence {
        const row = deduction.coset_in - 1;
        const col = deduction.gen;

        const item = &self.table_raw.items[row * self.num_gens + col];
        if (item.* != 0 and item.* != deduction.coset_out) {
            panic("Coincidence! C_{} == C_{}", .{ item.*, deduction.coset_out });
        }

        item.* = deduction.coset_out;
        return null;
    }

    pub fn is_full(self: Self) bool {
        for (self.table_raw.items) |item| {
            if (item == 0) {
                return false;
            }
        }
        return true;
    }

    pub fn find_unknown(self: Self) Unknown {
        for (0..self.num_cosets) |row_idx| {
            for (0..self.num_gens) |col_idx| {
                const coset_in = row_idx + 1;
                const gen = col_idx;

                if (self.table_raw.items[row_idx * self.num_gens + col_idx] == 0) {
                    return Unknown{
                        .gen = gen,
                        .coset_in = coset_in,
                    };
                }
            }
        }

        panic("If you're seeing this, :(", .{});
    }

    pub fn lookup(self: Self, coset_in: usize, gen: usize) ?usize {
        const row = coset_in - 1;
        const col = gen;

        if (coset_in > self.num_cosets) {
            panic("coset {} is too large (max {})", .{ coset_in, self.num_cosets });
        }
        if (gen >= self.num_gens) {
            panic("gen idx {} is too large (max {})", .{ gen, self.num_gens });
        }

        const coset_out = self.table_raw.items[row * self.num_gens + col];
        if (coset_out == 0) {
            return null;
        } else {
            return coset_out;
        }
    }
};

const RelTable = struct {
    const Self = @This();

    rel: Relation,
    num_cosets: usize,
    table_raw: ArrayList(usize),

    pub fn init(allocator: Allocator, rel: Relation) Self {
        return Self{
            .rel = rel,
            .num_cosets = 0,
            .table_raw = ArrayList(usize).init(allocator),
        };
    }

    pub fn deinit(self: Self) void {
        self.table_raw.deinit();
    }

    pub fn new_coset(self: *Self, id: usize) !void {
        for (0..self.rel.len - 1) |_| {
            try self.table_raw.append(0);
        }
        try self.table_raw.append(id);
        self.num_cosets += 1;
    }

    pub fn update_with_coset_table(self: *Self, coset_table: CosetTable, deductions: *ArrayList(Deduction)) !void {
        for (0..self.num_cosets) |row_idx| {
            var row = self.table_raw.items[row_idx * self.rel.len .. (row_idx + 1) * self.rel.len];
            for (0..self.rel.len - 1) |col_idx| {
                const gen = self.rel[col_idx];
                var coset_in = row_idx + 1;
                if (col_idx != 0) {
                    coset_in = row[col_idx - 1];
                }
                if (coset_in == 0) break;
                const coset_out = &row[col_idx];

                if (coset_out.* != 0) continue;

                if (coset_table.lookup(coset_in, gen)) |val| {
                    coset_out.* = val;

                    // adding the deduction completes the row iff this is the second last column
                    if (col_idx == self.rel.len - 2) {
                        // if so, we produce a new deduction
                        const new_deduction = Deduction{
                            .gen = self.rel[col_idx + 1],
                            .coset_in = val,
                            .coset_out = row[col_idx + 1],
                        };
                        try deductions.append(new_deduction);
                    }
                }
            }
        }
    }

    pub fn is_full(self: Self) bool {
        for (self.table_raw.items) |item| {
            if (item == 0) {
                return false;
            }
        }
        return true;
    }
};

/// Represents the fact that coset_in * gen = coset_out.
const Deduction = struct {
    gen: usize,
    coset_in: usize,
    coset_out: usize,
};

const Coincidence = struct {
    lower: usize,
    higher: usize,
};

const Unknown = struct {
    gen: usize,
    coset_in: usize,

    pub fn to_deduction(self: Unknown, coset_out: usize) Deduction {
        return Deduction{
            .coset_out = coset_out,
            .coset_in = self.coset_in,
            .gen = self.gen,
        };
    }
};

// Tests
const test_allocator = std.testing.allocator;
const expect = std.testing.expect;

test "symmetric group S3" {
    // The S3 group is presented by
    // <a, b | a^2, b^2, (ab)^3>
    var s3 = Presentation.init(test_allocator, 2);
    defer s3.deinit();
    try s3.add_relation(&[_]usize{ 0, 0 });
    try s3.add_relation(&[_]usize{ 1, 1 });
    try s3.add_relation(&[_]usize{ 0, 1, 0, 1, 0, 1 });

    const elements = try all_group_elements(test_allocator, s3);
    defer {
        for (elements) |ele| {
            test_allocator.free(ele);
        }
        test_allocator.free(elements);
    }

    std.debug.print("{any}\n", .{elements});

    try expect(elements.len == 6);

    try expect(std.mem.eql(usize, elements[0], &[_]usize{}));
    try expect(std.mem.eql(usize, elements[1], &[_]usize{0}));
    try expect(std.mem.eql(usize, elements[2], &[_]usize{1}));
    try expect(std.mem.eql(usize, elements[3], &[_]usize{ 1, 0, 1, 0 }));
    try expect(std.mem.eql(usize, elements[4], &[_]usize{ 1, 0 }));
    try expect(std.mem.eql(usize, elements[5], &[_]usize{ 1, 0, 1 }));
}

test "120 cell" {
    // This group is presented by
    // <a, b, c, d |
    //  a^2, b^2, c^2, d^2,
    //  (ac)^2, (ad)^2, (bd)^2,
    //  (ab)^5, (bc)^3, (cd)^3>
    var pres = Presentation.init(test_allocator, 4);
    defer pres.deinit();

    try pres.add_relation(&[_]usize{0} ** 2);
    try pres.add_relation(&[_]usize{1} ** 2);
    try pres.add_relation(&[_]usize{2} ** 2);
    try pres.add_relation(&[_]usize{3} ** 2);

    try pres.add_relation(&[_]usize{ 0, 2 } ** 2);
    try pres.add_relation(&[_]usize{ 0, 3 } ** 2);
    try pres.add_relation(&[_]usize{ 1, 3 } ** 2);

    try pres.add_relation(&[_]usize{ 0, 1 } ** 5);
    try pres.add_relation(&[_]usize{ 1, 2 } ** 3);
    try pres.add_relation(&[_]usize{ 2, 3 } ** 3);

    const elements = try all_group_elements(test_allocator, pres);
    defer {
        for (elements) |ele| {
            test_allocator.free(ele);
        }
        test_allocator.free(elements);
    }

    // std.debug.print("{any}\n", .{elements});

    try expect(elements.len == 600);
}
