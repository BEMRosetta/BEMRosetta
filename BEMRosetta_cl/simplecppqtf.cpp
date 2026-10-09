#include <Core/Core.h>

using namespace Upp;

namespace {

enum class TokenType {
    Normal,
    Keyword,
    String,
    Comment,
    Number,
    Preprocessor
};

bool IsCppIdentifierStart(int c) {
    return IsAlpha(c) || c == '_';
}

bool IsCppIdentifierChar(int c) {
    return IsAlNum(c) || c == '_';
}

bool IsCppKeyword(const String& word) {
    static const Index<String> keyword = [] {
        Index<String> k;

        const char *words[] = {
            "alignas", "alignof", "and", "and_eq", "asm", "auto", "bitand", "bitor", "bool", "break",
            "case", "catch", "char", "char8_t", "char16_t", "char32_t",
            "class", "compl", "concept", "const", "consteval",
            "constexpr", "constinit", "const_cast", "continue",
            "co_await", "co_return", "co_yield", "decltype", "default", "delete", "do", "double",
            "dynamic_cast", "else", "enum", "explicit", "export", 
            "extern", "false", "float", "for", "friend",
            "goto", "if", "inline", "int", "long", "mutable",
            "namespace", "new", "noexcept", "not", "not_eq",
            "nullptr", "operator", "or", "or_eq", "private", "protected", "public",
            "register", "reinterpret_cast", "requires", "return",
            "short", "signed", "sizeof", "static", "static_assert", "static_cast", "struct", "switch",
            "template", "this", "thread_local", "throw", "true", "try",
            "typedef", "typeid", "typename", "union", "unsigned", "using",
            "virtual", "void", "volatile", "wchar_t", "while", "xor", "xor_eq"
        };

        for (const char *word : words)
            k.Add(word);

        return k;
    }();

    return keyword.Find(word) >= 0;
}

String EscapeQtf(const String& text) {
    String out;

    for(int i = 0; i < text.GetCount(); ++i) {
        byte c = text[i];

        if(c == '\r')
            ;
        else if(c == '\n')
            out << '&';
        else if(c == '\t')
            out << "-|";
        else if(IsAlNum(c) ||
                c == ' ' ||
                c == '.' || c == ',' || c == ';' ||
                c == '!' || c == '?' || c == '%' ||
                c == '(' || c == ')' ||
                c == '/' || c == '<' || c == '>' ||
                c == '#' || c >= 128) {
            out.Cat(c);
        } else {
            out.Cat('`');
            out.Cat(c);
        }
    }
    return out;
}

String TokenColor(TokenType type) {
    switch(type) {
    case TokenType::Keyword:		return "(0.0.220)";       // Blue
    case TokenType::String:			return "(163.21.21)";     // Dark red
    case TokenType::Comment:		return "(0.128.0)";       // Green
    case TokenType::Number:			return "(128.0.128)";     // Purple
    case TokenType::Preprocessor:	return "(160.80.0)";      // Brown/orange
    default:						return Null;
    }
}

void AppendToken(String& qtf, const String& token, TokenType type) {
    if(token.IsEmpty())
        return;

    String escaped = EscapeQtf(token);

    if(type == TokenType::Normal)
        qtf << escaped;
    else
        qtf << "[@" << TokenColor(type) << " "  << escaped << "]";
}

bool IsNumberChar(int c) {
    return IsAlNum(c) || c == '.' || c == '\'' || c == '+' || c == '-';
}

}

String SimpleCppToQtf(const String& code) {
    String qtf;

    qtf << "[C ";

    int  i = 0;
    int  count = code.GetCount();
    bool beginning_of_line = true;

    while (i < count) {
        int start = i;
        int c = (byte)code[i];

        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            while (i < count) {
                int ch = (byte)code[i];

                if (ch != ' ' && ch != '\t' && ch != '\r' && ch != '\n')
                    break;

                if(ch == '\n')
                    beginning_of_line = true;

                ++i;
            }
            AppendToken(qtf, code.Mid(start, i - start),
                        TokenType::Normal);
            continue;
        }
        if (beginning_of_line && c == '#') {		// Preprocessor directive.
            ++i;

            while (i < count) {
                if (code[i] == '\n') {	// Continue if the previous character is '\'.
                    if (i > start && code[i - 1] == '\\') {
                        ++i;
                        continue;
                    }
                    break;
                }
                ++i;
            }
            AppendToken(qtf, code.Mid(start, i - start), TokenType::Preprocessor);

            beginning_of_line = false;
            continue;
        }
        beginning_of_line = false;

        if (c == '/' && i + 1 < count && code[i + 1] == '/') {		// Single-line comment
            i += 2;

            while (i < count && code[i] != '\n')
                ++i;

            AppendToken(qtf, code.Mid(start, i - start), TokenType::Comment);
            continue;
        }
        if (c == '/' && i + 1 < count && code[i + 1] == '*') {		// Block comment
            i += 2;

            while (i < count) {
                if (i + 1 < count && code[i] == '*' && code[i + 1] == '/') {
                    i += 2;
                    break;
                }

                if(code[i] == '\n')
                    beginning_of_line = true;

                ++i;
            }
            AppendToken(qtf, code.Mid(start, i - start), TokenType::Comment);
            continue;
        }
        if (c == '"' || c == '\'') {			// String or character literal
            char quote = c;
            ++i;

            while (i < count) {
                if (code[i] == '\\') {
                    i += min(2, count - i);
                    continue;
                }

                if (code[i] == quote) {
                    ++i;
                    break;
                }

                if (code[i] == '\n')
                    beginning_of_line = true;

                ++i;
            }
            AppendToken(qtf, code.Mid(start, i - start), TokenType::String);
            continue;
        }
        if (IsCppIdentifierStart(c)) {			// Identifier or keyword
            ++i;

            while(i < count &&
                  IsCppIdentifierChar((byte)code[i]))
                ++i;

            String word = code.Mid(start, i - start);

            AppendToken(qtf, word, IsCppKeyword(word) ? TokenType::Keyword : TokenType::Normal);
            continue;
        }
        if (IsDigit(c) || (c == '.' && i + 1 < count && IsDigit((byte)code[i + 1]))) {	// Numeric literal
            ++i;

            while (i < count &&
                  IsNumberChar((byte)code[i]))
                ++i;

            AppendToken(qtf, code.Mid(start, i - start), TokenType::Number);
            continue;
        }
        AppendToken(qtf, code.Mid(i, 1), TokenType::Normal);		// Operator or punctuation
        ++i;
    }
    qtf << "]";

    return qtf;
}

namespace {

bool IsPythonIdentifierStart(int c) {// Accept UTF-8 bytes without performing Unicode validation.
    return IsAlpha(c) || c == '_' || c >= 128;
}

bool IsPythonIdentifierChar(int c) {
    return IsPythonIdentifierStart(c) || IsDigit(c);
}

bool IsPythonKeyword(const String& word) {
    static const Index<String> keywords = [] {
        Index<String> k;

        const char *words[] = {
            "False", "None", "True",
            "and", "as", "assert", "async", "await",
            "break", "class", "continue", "def", "del",
            "elif", "else", "finally", "for",
            "global", "if", "in",
            "is", "lambda", "nonlocal", "not", "or",
            "pass", "raise", "return", "while",
            "with", "yield"
        };
        for (const char *word : words)
            k.Add(word);

        return k;
    }();

    return keywords.Find(word) >= 0;
}

bool IsPythonKeyword2(const String& word) {
    static const Index<String> keywords = [] {
        Index<String> k;

        const char *words[] = {
            "except", "from", "import", "try"
        };
        for (const char *word : words)
            k.Add(word);

        return k;
    }();

    return keywords.Find(word) >= 0;
}

bool IsPythonStringPrefix(const String& word) {
    String p = ToLower(word);

    return p == "r"  || p == "u"  || p == "b" ||
           p == "f"  || p == "t"  ||
           p == "br" || p == "rb" ||
           p == "fr" || p == "rf" ||
           p == "tr" || p == "rt";
}

// 'i' initially points to the opening quote, after any prefix.
void ScanPythonString(const String& code, int& i) {
    int count = code.GetCount();
    char quote = code[i];

    bool triple = i + 2 < count && code[i + 1] == quote && code[i + 2] == quote;

    i += triple ? 3 : 1;

    while (i < count) {
        if (code[i] == '\\') {		// Escaped quote, backslash, or physical newline.
            if(i + 2 < count && code[i + 1] == '\r' && code[i + 2] == '\n')
                i += 3;
            else
                i += min(2, count - i);

            continue;
        }
        if (code[i] == quote) {
            if(!triple) {
                ++i;
                return;
            }
            if(i + 2 < count && code[i + 1] == quote && code[i + 2] == quote) {
                i += 3;
                return;
            }
        }
        if (!triple && (code[i] == '\n' || code[i] == '\r'))	// Recover from an unterminated single-line string.
            return;

        ++i;
    }
}

bool IsPythonBaseDigit(int c, int base) {
    if (c >= '0' && c <= '9')
        return c - '0' < base;

    return base == 16 && ((c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F'));
}

void ScanPythonNumber(const String& code, int& i) {
    int count = code.GetCount();

    if (code[i] == '0' && i + 1 < count) {	// Binary, octal, or hexadecimal integer.
        int base = 0;

        switch (code[i + 1]) {
        case 'b': case 'B': base = 2;  break;
        case 'o': case 'O': base = 8;  break;
        case 'x': case 'X': base = 16; break;
        }
        if (base) {
            i += 2;

            while (i < count &&
                  (IsPythonBaseDigit((byte)code[i], base) ||
                   code[i] == '_'))
                ++i;
            return;
        }
    }
    while (i < count && (IsDigit((byte)code[i]) || code[i] == '_'))	// Decimal integer or floating-point mantissa.
        ++i;

    if (i < count && code[i] == '.') {
        ++i;
        while(i < count && (IsDigit((byte)code[i]) || code[i] == '_'))
            ++i;
    }
    if (i < count && (code[i] == 'e' || code[i] == 'E')) {			// Optional exponent. Signs belong here, not elsewhere.
        int exponent = i++;

        if(i < count && (code[i] == '+' || code[i] == '-'))
            ++i;

        if(i < count && IsDigit((byte)code[i])) {
            while(i < count &&
                  (IsDigit((byte)code[i]) || code[i] == '_'))
                ++i;
        }
        else
            i = exponent;
    }
    if (i < count && (code[i] == 'j' || code[i] == 'J'))		// Imaginary-number suffix.
        ++i;
}

}

String SimplePythonToQtf(const String& code) {
    String qtf;
    qtf << "[C ";

    int i = 0;
    int count = code.GetCount();

    while (i < count) {
        int start = i;
        int c = (byte)code[i];

        if(c == '#') {							// Comment.
            while(i < count && code[i] != '\n' && code[i] != '\r')
                ++i;

            AppendToken(qtf, code.Mid(start, i - start), TokenType::Comment);
            continue;
        }
        if (c == '\'' || c == '"') {				// String without a prefix.
            ScanPythonString(code, i);

            AppendToken(qtf, code.Mid(start, i - start),
                        TokenType::String);
            continue;
        }
        if (IsPythonIdentifierStart(c)) {			// Identifier, keyword, or string prefix.
            ++i;

            while (i < count && IsPythonIdentifierChar((byte)code[i]))
                ++i;

            String word = code.Mid(start, i - start);

            if (i < count && (code[i] == '\'' || code[i] == '"') && IsPythonStringPrefix(word)) {
                ScanPythonString(code, i);

                AppendToken(qtf, code.Mid(start, i - start),
                            TokenType::String);
            } else if (IsPythonKeyword(word))
                AppendToken(qtf, word, TokenType::Keyword);
            else if (IsPythonKeyword2(word))
                AppendToken(qtf, word, TokenType::Preprocessor);
            else 
                AppendToken(qtf, word, TokenType::Normal);

            continue;
        }
        if (IsDigit(c) || (c == '.' && i + 1 < count && IsDigit((byte)code[i + 1]))) {	// Numeric literal.
            ScanPythonNumber(code, i);

            AppendToken(qtf, code.Mid(start, i - start), TokenType::Number);
            continue;
        }
        AppendToken(qtf, code.Mid(i, 1), TokenType::Normal);	// Whitespace, operators, and punctuation.
        ++i;
    }

    qtf << "]";
    return qtf;
}