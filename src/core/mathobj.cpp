#include "parser.h"
#include "mathobj.h"

//SW: This calls the function with the opcode after it's been parsed
_PMathObj _MathObject::Execute (long opCode, _PMathObj p, _PMathObj p2)   // execute this operation with the second arg if necessary
{
    switch (opCode) {
    case HY_OP_CODE_NOT: // !
        return LNot();
        break;
    case HY_OP_CODE_NEQ: // !=
        return NotEqual(p);
        break;
    case HY_OP_CODE_IDIV: // $
        return longDiv(p);
        break;
    case HY_OP_CODE_MOD: // %
        return lDiv(p);
        break;
    case HY_OP_CODE_AND: // &&
        return LAnd(p);
        break;
    case HY_OP_CODE_MUL: // *
        return Mult(p);
        break;
    case HY_OP_CODE_ADD: // +
        if (p) {
            return Add(p);
        } else {
            return Sum ();
        }
        break;
    case HY_OP_CODE_SUB: // -
        if (p) {
            return Sub(p);
        } else {
            return Minus();
        }
        break;
    case HY_OP_CODE_DIV: // /
        return Div(p);
        break;
    case HY_OP_CODE_LESS: // <
        return Less(p);
        break;
    case HY_OP_CODE_LEQ: // <=
        return LessEq(p);
        break;
    case HY_OP_CODE_EQ: // ==
        return AreEqual(p);
        break;
    case HY_OP_CODE_GREATER: // >
        return Greater(p);
        break;
    case HY_OP_CODE_GEQ: // >=
        return GreaterEq(p);
        break;
    case HY_OP_CODE_ABS: // Abs
        return Abs();
        break;
    case HY_OP_CODE_ARCTAN: // Arctan
        return Arctan();
        break;
    case HY_OP_CODE_BETA: // Beta
        return Beta(p);
        break;
    case HY_OP_CODE_CCHI2: // CChi2
        return CChi2(p);
        break;
    case HY_OP_CODE_CGAMMADIST: // CGammaDist
        return CGammaDist(p,p2);
        break;
    case HY_OP_CODE_COS: // Cos
        return Cos();
        break;
    case HY_OP_CODE_COLUMNS: // Columns
    case HY_OP_CODE_ROWS: // Rows
        return new _Constant (0.0);
        break;
    case HY_OP_CODE_ERF: // Erf
        return Erf();
        break;
    case HY_OP_CODE_EXP: // Exp
        return Exp();
        break;
    case HY_OP_CODE_FORMAT: // Format
        return FormatNumberString(p,p2);
        break;
    case HY_OP_CODE_GAMMA: // Gamma
        return Gamma();
        break;
    case HY_OP_CODE_GAMMADIST: // GammaDist
        return GammaDist(p,p2);
        break;
    case HY_OP_CODE_IBETA: // IBeta
        return IBeta(p,p2);
        break;
    case HY_OP_CODE_IGAMMA: // IGamma
        return IGamma(p);
        break;
    case HY_OP_CODE_INVCHI2: // InvChi2
        return InvChi2(p);
        break;
    case HY_OP_CODE_LNGAMMA: // LnGamma
        return LnGamma();
        break;
    case HY_OP_CODE_LOG: // Log
        return Log();
        break;
    case HY_OP_CODE_MAX: // Max
        return Max(p);
        break;
    case HY_OP_CODE_MIN: // Min
        return Min(p);
        break;
    case HY_OP_CODE_RANDOM: // Random
        return Random(p);
        break;
    case HY_OP_CODE_SIMPLEX: // Simplex
        return Simplex();
        break;
    case HY_OP_CODE_SIN: // Sin
        return Sin();
        break;
    case HY_OP_CODE_SQRT: // Sqrt
        return Sqrt();
        break;
    case HY_OP_CODE_TAN: // Tan
        return Tan();
        break;
    case HY_OP_CODE_TIME: // Time
        return Time();
        break;
    case HY_OP_CODE_TYPE: // Type
        return Type();
        break;
    case HY_OP_CODE_ZCDF: // ZCDF
        return ZCDF();
        break;
    case HY_OP_CODE_POWER: // ^
        return Raise(p);
        break;
    case HY_OP_CODE_OR: // ||
        return LOr(p);
        break;
    default: {
        WarnNotDefined (this, opCode);
        return new _Constant (0.0);
    }
    }
    return (_PMathObj)makeDynamic();
}

//__________________________________________________________________________________

BaseRef _MathObject::makeDynamic (void)
{
    return(_PMathObj)checkPointer(new _MathObject);
}

//SW: Why do we need a string for the type?
_PMathObj _MathObject::Type (void)
{
    _FString * ts = new _FString();
    switch (ObjectClass()) {

    case NUMBER:
        *(ts->theString)="Number";
        break;
    case MATRIX:
        *(ts->theString)="Matrix";
        break;
    case CONTAINER:
        *(ts->theString)="Container";
        break;
    case TREE_NODE:
        *(ts->theString)="TreeNode";
        break;
    case TREE:
        *(ts->theString)="Tree";
        break;
    case STRING:
        *(ts->theString)="String";
        break;
    case ASSOCIATIVE_LIST:
        *(ts->theString)="AssociativeList";
        break;
    case TOPOLOGY:
        *(ts->theString)="Topology";
        break;
    case POLYNOMIAL:
        *(ts->theString)="Polynomial";
        break;
    default:
        *(ts->theString) = "Unknown";

    }

    return ts;
}
