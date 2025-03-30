//! Implementation of Todd-Coxeter algorithm to find unique elements of a
//! group given by a finite presentation (`all_group_elements`).

const std = @import("std");
const Allocator = std.mem.Allocator;
const ArrayList = std.ArrayList;
const AutoHashMap = std.AutoHashMap;
const panic = std.debug.panic;
const log = std.log;
const test_allocator = std.testing.allocator;
const expect = std.testing.expect;

const presentation = @import("presentation.zig");
const Presentation = presentation.Presentation;
const Relation = presentation.Relation;

/// Find all the distinct elements of the group represented
/// by the given presentation. All generators must be involutions.
/// The return elements are represented by a word made from the presentation's generators.
pub fn all_group_elements(allocator: Allocator, pres: Presentation) !Elements {
    var coset_table = CosetTable.init(allocator, pres.num_gens);
    defer coset_table.deinit();

    var rel_tables = ArrayList(RelTable).init(allocator);
    defer {
        for (rel_tables.items) |*rel_table| {
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

    // Deductions and coincidences which have yet to be procesed are stored in a stack
    var deductions = ArrayList(Deduction).init(allocator);
    defer deductions.deinit();
    var coincidences = ArrayList(Coincidence).init(allocator);
    defer coincidences.deinit();

    // TODO: special case if any of the rel_tables have a width (relation length) of 1.
    // The row is already completed and so we have a new deduction.

    var num_cosets: usize = 1;
    while (true) {
        std.log.debug("Coset table:\n{}\n", .{coset_table});
        // handle any deductions we have
        while (deductions.pop()) |deduction| {
            std.log.debug("Adding deduction {any}\n", .{deduction});
            if (coset_table.add_deduction(deduction)) |coincidence| {
                try coincidences.append(coincidence);
            }

            while (coincidences.pop()) |coincidence| {
                // If filling the coset table produces a coincidence Ci = Cj, i < j
                // we need to replace all the instances of j with i.
                // Then fill in all possible entries of the tables

                std.log.debug("Coincidence! {}\n", .{coincidence});
                try coset_table.handle_coincidence(coincidence, &coincidences);
                for (rel_tables.items) |*rel_table| {
                    try rel_table.handle_coincidence(coincidence, &deductions);
                }
                for (deductions.items) |*d| {
                    if (d.coset_in == coincidence.higher) {
                        d.coset_in = coincidence.lower;
                    }
                    if (d.coset_out == coincidence.higher) {
                        d.coset_out = coincidence.lower;
                    }
                }
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
        std.log.debug("Adding coset {}\n", .{num_cosets});
        const unknown = coset_table.find_unknown();
        try coset_table.new_coset(num_cosets);
        for (rel_tables.items) |*rel_table| {
            try rel_table.new_coset(num_cosets);
        }
        try deductions.append(unknown.to_deduction(num_cosets));
    }

    std.log.debug("Coset table:\n{}\n", .{coset_table});

    // Now that the coset table is complete, we can build the group's elements
    return try elements_from_coset_table(allocator, coset_table);
}

fn elements_from_coset_table(allocator: Allocator, coset_table: CosetTable) !Elements {
    var elements = try allocator.alloc(ArrayList(usize), coset_table.coset_count());
    for (0..coset_table.coset_count()) |i| {
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
                const word = elements[coset_in - 1];
                try elements[coset_out - 1].appendSlice(word.items);
                try elements[coset_out - 1].append(gen);
                try coset_stack.append(coset_out);
            }
        }
    }

    var non_opt_elements = try allocator.alloc(GeneratorWord, coset_table.coset_count());
    for (0..coset_table.coset_count()) |i| {
        non_opt_elements[i] = try elements[i].toOwnedSlice();
    }

    return Elements{
        .items = non_opt_elements,
        .allocator = allocator,
    };
}

pub const GeneratorWord = []usize;

const Elements = struct {
    items: []GeneratorWord,
    allocator: Allocator,

    pub fn deinit(self: Elements) void {
        for (self.items) |word| {
            self.allocator.free(word);
        }
        self.allocator.free(self.items);
    }
};

const CosetTable = struct {
    const Self = @This();

    const Key = struct { coset: usize, gen: usize };
    const Map = AutoHashMap(Key, ?usize);

    num_gens: usize,
    map: Map,
    /// This gives our unknowns an order for `find_unknown` to use
    next_unknowns: ArrayList(Key),
    allocator: Allocator,

    pub fn init(allocator: Allocator, num_gens: usize) Self {
        return Self{
            .num_gens = num_gens,
            .map = Map.init(allocator),
            .next_unknowns = ArrayList(Key).init(allocator),
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        self.map.deinit();
        self.next_unknowns.deinit();
    }

    pub fn coset_count(self: Self) usize {
        return self.map.count() / self.num_gens;
    }

    pub fn new_coset(self: *Self, id: usize) !void {
        std.log.info("new coset: {}", .{id});
        for (0..self.num_gens) |gen| {
            const key: Key = .{ .coset = id, .gen = gen };
            try self.map.put(key, null);
            try self.next_unknowns.append(key);
        }
    }

    pub fn add_deduction(self: *Self, deduction: Deduction) ?Coincidence {
        const key: Key = .{ .coset = deduction.coset_in, .gen = deduction.gen };
        const item_ptr = self.map.getPtr(key).?;
        if (item_ptr.*) |item| {
            if (item != deduction.coset_out) {
                if (item < deduction.coset_out) {
                    return Coincidence{
                        .lower = item,
                        .higher = deduction.coset_out,
                    };
                } else {
                    return Coincidence{
                        .higher = item,
                        .lower = deduction.coset_out,
                    };
                }
            }
        }

        item_ptr.* = deduction.coset_out;
        return null;
    }

    pub fn is_full(self: Self) bool {
        var iter = self.map.valueIterator();
        while (iter.next()) |value| {
            if (value.* == null) {
                return false;
            }
        }
        return true;
    }

    fn is_ignored(self: Self, coset: usize) bool {
        return std.mem.indexOfScalar(usize, self.ignored_rows.items, coset) != null;
    }

    pub fn find_unknown(self: *Self) Unknown {
        // Search `next_unkowns` for first actual unknown
        while (true) {
            // TODO: this remove should be O(1)
            const key = self.next_unknowns.orderedRemove(0);
            if (self.lookup(key.coset, key.gen) == null) {
                return Unknown{
                    .coset_in = key.coset,
                    .gen = key.gen,
                };
            }
        }
    }

    pub fn lookup(self: Self, coset_in: usize, gen: usize) ?usize {
        const key: Key = .{ .coset = coset_in, .gen = gen };
        return self.map.get(key) orelse {
            panic("key {} {} not found", .{ coset_in, gen });
        };
    }

    pub fn handle_coincidence(self: *Self, coin: Coincidence, coincidences: *ArrayList(Coincidence)) !void {
        const Entry = struct { key: Key, value: ?usize };
        var to_update = ArrayList(Entry).init(self.allocator);
        defer to_update.deinit();

        // Find entries that need to be updated
        var iter = self.map.iterator();
        while (iter.next()) |entry| {
            if (entry.key_ptr.*.coset == coin.higher or entry.value_ptr.* == coin.higher) {
                try to_update.append(.{ .key = entry.key_ptr.*, .value = entry.value_ptr.* });
            }
        }

        for (to_update.items) |*entry| {
            const old_opt = self.map.get(entry.key).?;
            _ = self.map.remove(entry.key);
            if (old_opt != null and entry.value != null and old_opt.? != entry.value.?) {
                // If the old value existed, and the new value exists, and both are different
                const new_coin = Coincidence{
                    .lower = @min(old_opt.?, entry.value.?),
                    .higher = @max(old_opt.?, entry.value.?),
                };
                try coincidences.append(new_coin);
            }
            if (entry.key.coset == coin.higher) {
                entry.key.coset = coin.lower;
            }
            if (entry.value == coin.higher) {
                entry.value = coin.lower;
            }
            try self.map.put(entry.key, entry.value);
        }
    }

    pub fn format(self: Self, comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
        _ = fmt;
        _ = options;

        try writer.print("\t", .{});
        for (0..self.num_gens) |gen| {
            try writer.print("g{}\t", .{gen});
        }
        try writer.print("\n", .{});

        var processed_cosets: usize = 0;
        var coset_id: usize = 0;

        while (processed_cosets * self.num_gens < self.map.count()) {
            if (!self.map.contains(.{ .coset = coset_id, .gen = 0 })) {
                coset_id += 1;
                continue;
            }

            try writer.print("{}\t", .{coset_id});
            for (0..self.num_gens) |gen| {
                if (self.map.get(.{ .coset = coset_id, .gen = gen })) |coset_out_opt| {
                    if (coset_out_opt) |coset_out| {
                        try writer.print("{}\t", .{coset_out});
                    } else {
                        try writer.print("?\t", .{});
                    }
                } else {
                    panic("bad state", .{});
                }
            }
            try writer.print("\n", .{});
            coset_id += 1;
            processed_cosets += 1;
        }
    }
};

const RelTable = struct {
    const Self = @This();

    const Key = struct { coset: usize };
    const Table = ArrayList([]?usize);

    rel: Relation,
    table: Table,
    allocator: Allocator,

    pub fn init(allocator: Allocator, rel: Relation) Self {
        return Self{
            .rel = rel,
            .table = Table.init(allocator),
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        for (self.table.items) |row| {
            self.allocator.free(row);
        }
        self.table.deinit();
    }

    pub fn new_coset(self: *Self, id: usize) !void {
        var row = try self.allocator.alloc(?usize, self.rel.len + 1);
        errdefer self.allocator.free(row);

        row[0] = id;
        for (1..self.rel.len) |i| {
            row[i] = null;
        }
        row[self.rel.len] = id;

        try self.table.append(row);
    }

    pub fn update_with_coset_table(self: *Self, coset_table: CosetTable, deductions: *ArrayList(Deduction)) !void {
        var i: usize = 0;
        while (i < self.table.items.len) {
            const row = self.table.items[i];
            if (try self.update_row_with_coset_table(row, coset_table, deductions)) {
                // since the row was completed, we can remove it
                _ = self.table.orderedRemove(i);
            } else {
                i += 1;
            }
        }
    }

    /// Returns a boolean indicated whether the row was completed
    fn update_row_with_coset_table(self: *Self, row: []?usize, coset_table: CosetTable, deductions: *ArrayList(Deduction)) !bool {
        // Do two passes, once LTR and once RTL
        for (0..self.rel.len) |col_idx| {
            const gen = self.rel[col_idx];
            const left = &row[col_idx];
            const right = &row[col_idx + 1];

            if (left.* != null and right.* == null) {
                // can potentially update right value
                if (coset_table.lookup(left.*.?, gen)) |val| {
                    right.* = val;
                    // adding the deduction completes the row if the value *after* right is not null
                    if (row[col_idx + 2] != null) {
                        // if so, we produce a new deduction
                        const new_deduction = Deduction{
                            .gen = self.rel[col_idx + 1],
                            .coset_in = val,
                            .coset_out = row[col_idx + 2].?,
                        };
                        try deductions.append(new_deduction);
                        return true;
                    }
                }
            }
        }

        for (0..self.rel.len) |col_idx_rev| {
            const col_idx = self.rel.len - 1 - col_idx_rev;
            const gen = self.rel[col_idx];
            const left = &row[col_idx];
            const right = &row[col_idx + 1];

            if (left.* == null and right.* != null) {
                // can potentially update left value
                if (coset_table.lookup(right.*.?, gen)) |val| {
                    left.* = val;
                    // adding the deduction completes the row if the value *before* left is not null
                    if (row[col_idx - 1] != null) {
                        // if so, we produce a new deduction
                        const new_deduction = Deduction{
                            .gen = self.rel[col_idx - 1],
                            .coset_in = val,
                            .coset_out = row[col_idx - 1].?,
                        };
                        try deductions.append(new_deduction);
                        return true;
                    }
                }
            }
        }

        return false;
    }

    pub fn is_full(self: Self) bool {
        for (self.table.items) |row| {
            for (row) |e| {
                if (e == null) return false;
            }
        }
        return true;
    }

    pub fn handle_coincidence(self: *Self, coin: Coincidence, deductions: *ArrayList(Deduction)) !void {
        _ = deductions;

        for (self.table.items) |row| {
            for (row) |*coset_out| {
                if (coset_out.* == coin.higher) {
                    coset_out.* = coin.lower;
                }
            }
        }
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
fn elements_contains(elements: Elements, word: []const usize) bool {
    for (elements.items) |element| {
        if (std.mem.eql(usize, element, word)) {
            return true;
        }
    }

    return false;
}

test "symmetric group S3" {
    // The S3 group is presented by
    // <a, b | a^2, b^2, (ab)^3>
    var s3 = Presentation.init(test_allocator, 2);
    defer s3.deinit();
    try s3.add_relation(&[_]usize{ 0, 0 });
    try s3.add_relation(&[_]usize{ 1, 1 });
    try s3.add_relation(&[_]usize{ 0, 1, 0, 1, 0, 1 });

    const elements = try all_group_elements(test_allocator, s3);
    defer elements.deinit();

    std.debug.print("{any}\n", .{elements});

    try expect(elements.items.len == 6);

    try expect(elements_contains(elements, &[_]usize{}));
    try expect(elements_contains(elements, &[_]usize{0}));
    try expect(elements_contains(elements, &[_]usize{1}));
    try expect(elements_contains(elements, &[_]usize{ 1, 0, 1, 0 }));
    try expect(elements_contains(elements, &[_]usize{ 1, 0 }));
    try expect(elements_contains(elements, &[_]usize{ 1, 0, 1 }));
}

test "square" {
    // This group is presented by
    // <a, b | a^2, b^2, (ab)^4>
    var pres = Presentation.init(test_allocator, 2);
    defer pres.deinit();

    try pres.add_relation(&[_]usize{0} ** 2);
    try pres.add_relation(&[_]usize{1} ** 2);

    try pres.add_relation(&[_]usize{ 0, 1 } ** 4);

    const elements = try all_group_elements(test_allocator, pres);
    defer elements.deinit();

    std.debug.print("{any}\n", .{elements});

    try expect(elements.items.len == 8);
}

test "cube" {
    // This group is presented by
    // <a, b, c |
    //  a^2, b^2, c^2, (ac)^2, (ab)^4, (bc)^3>
    var pres = Presentation.init(test_allocator, 3);
    defer pres.deinit();

    try pres.add_relation(&[_]usize{0} ** 2);
    try pres.add_relation(&[_]usize{1} ** 2);
    try pres.add_relation(&[_]usize{2} ** 2);
    try pres.add_relation(&[_]usize{ 0, 2 } ** 2);
    try pres.add_relation(&[_]usize{ 0, 1 } ** 4);
    try pres.add_relation(&[_]usize{ 1, 2 } ** 3);

    const elements = try all_group_elements(test_allocator, pres);
    defer elements.deinit();

    try expect(elements.items.len == 48);
}

test "[5, 3, 3]" {
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
    defer elements.deinit();

    try expect(elements.items.len == 14400); // absolutely wild number ngl
}
