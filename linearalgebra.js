Numbas.addExtension('linearalgebra',['jme','jme-display'],function(extension) {

var math = Numbas.math;

var Fraction = extension.Fraction = function(n) {
    if(typeof(n)=='number') {
        var c = math.rationalApproximation(n);
        this.n = c[0];
        this.d = c[1];
    } else {
        this.n = n.n;
        this.d = n.d;
    }
    this.tidy();
}
Fraction.prototype = {
    toString: function() {
        return this.d==1 ? this.n+'' : this.n+'/'+this.d;
    },

    toLaTeX: function() {
        return this.d==1 ? this.n+'' : (this.n<0 ? '-': '')+'\\frac{'+(Math.abs(this.n))+'}{'+this.d+'}';
    },

    /** Ensure fraction is reduced, and denominator is positive
     */
    tidy: function() {
        if(this.n==0) {
            this.n = 0;
            this.d = 1;
        }
        var g = math.gcd(this.n,this.d) * (this.d<0 ? -1 : 1);
        this.n /= g;
        this.d /= g;
    },

    add: function(f2) {
        return new Fraction({n: this.n*f2.d + f2.n*this.d, d: this.d*f2.d});
    },
    sub: function(f2) {
        return new Fraction({n: this.n*f2.d - f2.n*this.d, d: this.d*f2.d});
    },
    mul: function(f2) {
        return new Fraction({n: this.n*f2.n, d: this.d*f2.d});
    },
    div: function(f2) {
        return new Fraction({n: this.n*f2.d, d: this.d*f2.n});
    },
    is_zero: function() {
        return this.n==0;
    },
    is_one: function() {
        return this.n==this.d;
    },
    reciprocal: function() {
        return new Fraction({n:this.d,d:this.n});
    }
}

var fraction_matrix = extension.fraction_matrix = function(matrix) {
    var o = matrix.map(function(r){return r.map(function(c){ return new Fraction(c)})});
    o.rows = matrix.rows;
    o.columns = matrix.columns;
    return o;
}
var unfraction_matrix = extension.unfraction_matrix = function(matrix) {
    var o = matrix.map(function(r){return r.map(function(c){return c.n/c.d})});
    o.rows = matrix.rows;
    o.columns = matrix.columns;
    return o;
}

function logger(operations,matrix) {
    return function log(message,include_matrix,options) {
        include_matrix = include_matrix===undefined ? true : include_matrix;
        var lmatrix;
        if(include_matrix) {
            lmatrix = matrix.map(function(r){return r.slice()});
            lmatrix.rows = matrix.rows;
            lmatrix.columns = matrix.columns;
        }
        var l = {message:message, matrix: lmatrix};
        for(var key in options) {
            l[key] = options[key];
        }
        operations.push(l);
    }
}

var row_echelon_form = function(matrix) {
    /** Put a matrix representing a system of equations in row-echelon form.
    * Can:
    * * Swap two rows
    * * Multiply a row by a scalar
    * * Subtract a multiple of one row from another
    * For each row of the output, the first non-zero entry is 1, and strictly to the right of the first non-zero entry in the row above.
    * Works over the rationals: input is a matrix of objects {n: numerator,d: denominator}.
    * Output is an object {matrix, operations}, where operations is a list of descriptions of each step of the process, of the form {message: string, matrix: the state of the matrix after the operation}.
    */
    matrix = matrix.map(function(r){return r.slice()});
    var rows = matrix.length;
    var columns = matrix[0].length;
    matrix.rows = rows;
    matrix.columns = columns;
    
    var operations = [];
    var log = logger(operations,matrix);

    var current_row = 0;
    // for each column, there should be at most one row with a 1 in that column, and every other row should have 0 in that column
    for(var leader_column=0;leader_column<columns;leader_column++) {
        // find the first row with a non-zero in that column
        for(var row=current_row;row<rows;row++) {
            if(!matrix[row][leader_column].is_zero()) {
                break;
            }
        }
        // if we found a row with a non-zero in the leader column 
        if(row<rows) {
            // swap that row with the <current_row>th one
            if(row!=current_row) {
                var tmp = matrix[row];
                matrix[row] = matrix[current_row];
                matrix[current_row] = tmp;
                log("Row "+(row+1)+" has a non-zero entry in column "+(leader_column+1)+", so it should go before row "+(current_row+1)+". Swap row "+(row+1)+" with row "+(current_row+1)+".",true,{determinant_scale:new Fraction(-1)});
            }

            // multiply this row so the leader column has a 1 in it
            var leader = matrix[current_row][leader_column];
            if(!leader.is_one()) {
                matrix[current_row] = matrix[current_row].map(function(c){ return c.div(leader)});
                log("Divide row "+(current_row+1)+" by \\("+leader+"\\), so that the first non-zero entry is \\(1\\).",true, {determinant_scale:leader.reciprocal()});
            }

            // subtract multiples of this row from every other row so they all have a zero in this column
            var sub = function(a,b){ return a.sub(b); };
            var add = function(a,b){ return a.add(b); };
            for(var row=current_row+1;row<rows;row++) {
                if(row!=current_row && !matrix[row][leader_column].is_zero()) {
                    var scale = matrix[row][leader_column];
                    var op = sub;
                    if(scale.n<0) {
                        scale = new Fraction({n:-scale.n,d:scale.d});
                        op = add;
                    }
                    matrix[row] = matrix[row].map(function(c,i) { 
                        var res = op(c,matrix[current_row][i].mul(scale));
                        return res;
                    });
                    var mop = op==sub ? "Subtract" : "Add";
                    var mverb = op==sub ? "from" : "to";
                    log(mop+" "+(scale.is_one() ? "" : "\\("+scale+"\\) times ")+"row "+(current_row+1)+" "+mverb+" row "+(row+1)+".");
                }
            }
            current_row += 1;
        }
    }
    if(operations.length>0) {
        log("The matrix is now in row echelon form.",false);
    }
    return {
        matrix: matrix,
        operations: operations
    };
}
extension.row_echelon_form = function(matrix) {
    matrix = fraction_matrix(matrix);
    var res = row_echelon_form(matrix);
    res.matrix = unfraction_matrix(res.matrix);
    return res;
}

var reduced_row_echelon_form = function(matrix) {
    /** Put a matrix representing a system of equations in reduced row-echelon form.
     * Can:
     * * Swap two rows
     * * Multiply a row by a scalar
     * * Subtract a multiple of one row from another
     * As well as being in row-echelon form, the matrix has the property that the first non-zero entry in each row is also the only non-zero entry in its column.
     * Works over the rationals: input is a matrix of objects {n: numerator,d: denominator}.
     * Output is an object {matrix, operations}, where operations is a list of descriptions of each step of the process, of the form {message: string, matrix: the state of the matrix after the operation}.
     */
    matrix = matrix.map(function(r){return r.slice()});
    var res = row_echelon_form(matrix);
    matrix = res.matrix;
    var operations = res.operations.slice();

    var rows = matrix.length;
    var columns = matrix[0].length;
    matrix.rows = rows;
    matrix.columns = columns;

    var log = logger(operations,matrix);

    var sub = function(a,b){ return a.sub(b); };
    var add = function(a,b){ return a.add(b); };

    for(var row=0;row<rows;row++) {
        for(var column=0;column<columns && matrix[row][column].is_zero();column++) {}
        
        if(column==columns) {
            continue;
        }
        for(var vrow = 0;vrow<rows;vrow++) {
            if(vrow!=row && !matrix[vrow][column].is_zero()) {
                
                var scale = matrix[vrow][column];
                if(!scale.is_zero()) {
                    var op = sub;
                    if(scale.n<0) {
                        op = add;
                        scale = new Fraction({n:-scale.n, d:scale.d});
                    }
                    matrix[vrow] = matrix[vrow].map(function(c,i) { 
                        return op(c,matrix[row][i].mul(scale));
                    });

                    var mop = op==sub ? "subtract" : "add";
                    var mverb = op==sub ? "from" : "to";
                    log("We want a zero in column "+(column+1)+" of row "+(vrow+1)+": "+mop+" "+(scale.is_one() ? "" : "\\("+scale+"\\) times ")+"row "+(row+1)+" "+mverb+" row "+(vrow+1)+".");
                }
            }
        }
    }
    if(operations.length>0) {
        log("The matrix is now in reduced row echelon form.",false);
    }
    return {
        matrix: matrix,
        operations: operations
    };
}
extension.reduced_row_echelon_form = function(matrix) {
    matrix = fraction_matrix(matrix);
    var res = reduced_row_echelon_form(matrix);
    res.matrix = unfraction_matrix(res.matrix);
    return res;
}

/** Is the given matrix in row echelon form?
 * If not, throws an error with an explanation why it isn't.
 */
var is_row_echelon_form = extension.is_row_echelon_form = function(matrix) {
    var leader = -1;
    var rows = matrix.length;
    var columns = matrix[0].length;
    for(var row=0;row<rows;row++) {
        for(var column=0;column<columns;column++) {
            var cell = matrix[row][column];
            if(column<=leader) {
                if(!cell.is_zero()) {
                    throw(new Error("The first non-zero entry in row "+(row+1)+" is not strictly to the right of the first non-zero entries in the rows above."));
                } 
            } else {
                leader = column;
            }
        }
    }
    return true;
}
extension.is_row_echelon_form = function(matrix) {
    matrix = fraction_matrix(matrix);
    return is_row_echelon_form(matrix);
}

/** Is the given matrix in row echelon form?
 * If not, throws an error with an explanation why it isn't.
 */
var is_reduced_row_echelon_form = extension.is_reduced_row_echelon_form = function(matrix) {
    is_row_echelon_form(matrix); // this will throw an error if the matrix is not in row echelon form

    for(var row=0;row<matrix.rows;row++) {
        for(var column=0;column<matrix.columns;column++) {
            var cell = matrix[row][column];
            if(!cell.is_zero()) {
                if(!cell.is_one()) {
                    throw(new Error("The first non-zero entry in row "+(row+1)+" is not 1."))
                }
                for(var vrow=0;vrow<matrix.rows;vrow++) {
                    if(vrow!=row && !matrix[vrow][column].is_zero()) {
                        throw(new Error("There is more than one non-zero value in column "+(column+1)+"."));
                    }
                }
                break;
            }
        }
    }
    return true;
}
extension.is_reduced_row_echelon_form = function(matrix) {
    matrix = fraction_matrix(matrix);
    return is_reduced_row_echelon_form(matrix);
}

var scope = extension.scope;
var jme = Numbas.jme;
var funcObj = jme.funcObj;
var TMatrix = jme.types.TMatrix;
var TString = jme.types.TString;
var THTML = jme.types.THTML;
var TBool = jme.types.TBool;

function element(name,attrs,content) {
    var e = document.createElement(name);
    for(var k in attrs) {
        e.setAttribute(k,attrs[k]);
    }
    if(content!==undefined) {
        e.innerHTML = content;
    }
    return e;
}

scope.addFunction(new funcObj('row_echelon_form',[TMatrix],TMatrix,function(matrix) {
    matrix = fraction_matrix(matrix);
    var res = row_echelon_form(matrix);
    var omatrix = unfraction_matrix(res.matrix);
    return omatrix;
}));

function show_steps(steps,describe_determinant) {
    var ops = element('ul');
    if(describe_determinant) {
        var li = element('li',{},'Let the determinant of the matrix be \\(d\\)');
    }
    var d = new Fraction(1);
    steps.map(function(o) {
        var li = element('li');
        li.appendChild(element('span',{},o.message));
        if(o.matrix) {
            var m = new TMatrix(unfraction_matrix(o.matrix));
            li.appendChild(element('span',{},'\\['+jme.display.texify({tok:m},{fractionnumbers:true})+'\\]'));
            if(describe_determinant && o.determinant_scale) {
                d = d.mul(o.determinant_scale);
                li.appendChild(element('p',{},'The determinant of this matrix is \\('+(Math.abs(d.n)==d.d ? d.n<0 ? '-' : '' : d.toLaTeX())+' d\\).'));
            }
        }
        ops.appendChild(li);
    });
    return new THTML(ops);
}

scope.addFunction(new funcObj('row_echelon_form_display',[TMatrix],THTML,function(matrix) {
    matrix = fraction_matrix(matrix);
    var res = row_echelon_form(matrix);
    return show_steps(res.operations);
},{unwrapValues:true}));

scope.addFunction(new funcObj('row_echelon_form_display_determinant',[TMatrix],THTML,function(matrix) {
    matrix = fraction_matrix(matrix);
    var res = row_echelon_form(matrix);
    return show_steps(res.operations,true);
},{unwrapValues:true}));

scope.addFunction(new funcObj('is_row_echelon_form',[TMatrix],TBool,function(matrix) {
    matrix = fraction_matrix(matrix);
    try {
        return is_row_echelon_form(matrix);
    } catch(e) {
        return false;
    }
}));

scope.addFunction(new funcObj('describe_why_row_echelon_form',[TMatrix],TString,function(matrix) {
    matrix = fraction_matrix(matrix);
    try {
        is_row_echelon_form(matrix);
        return "The matrix is in row echelon form.";
    } catch(e) {
        return e.message;
    }
}));

scope.addFunction(new funcObj('reduced_row_echelon_form',[TMatrix],TMatrix,function(matrix) {
    matrix = fraction_matrix(matrix);
    var res = reduced_row_echelon_form(matrix);
    var omatrix = unfraction_matrix(res.matrix);
    return omatrix;
}));

scope.addFunction(new funcObj('reduced_row_echelon_form_display',[TMatrix],THTML,function(matrix) {
    matrix = fraction_matrix(matrix);
    var res = reduced_row_echelon_form(matrix);
    return show_steps(res.operations);
},{unwrapValues:true}));

scope.addFunction(new funcObj('is_reduced_row_echelon_form',[TMatrix],TBool,function(matrix) {
    matrix = fraction_matrix(matrix);
    try {
        return is_reduced_row_echelon_form(matrix);
    } catch(e) {
        return false;
    }
}));

scope.addFunction(new funcObj('describe_why_reduced_row_echelon_form',[TMatrix],TString,function(matrix) {
    matrix = fraction_matrix(matrix);
    try {
        is_reduced_row_echelon_form(matrix);
        return "The matrix is in reduced row echelon form.";
    } catch(e) {
        return e.message;
    }
}));

});
