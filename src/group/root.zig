const presentation = @import("presentation.zig");
pub const Presentation = presentation.Presentation;
pub const Relation = presentation.Relation;
const todd_coxeter = @import("todd-coxeter.zig");
pub const all_group_elements = todd_coxeter.all_group_elements;

test {
    _ = presentation;
    _ = todd_coxeter;
}
