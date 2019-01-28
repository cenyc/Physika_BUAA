/*
 * @file linear_system_solver.h
 * @brief base class of all solvers to solve linear system Ax = b
 * @author Fei Zhu
 *
 * This file is part of Physika, a versatile physics simulation library.
 * Copyright (C) 2013- Physika Group.
 *
 * This Source Code Form is subject to the terms of the GNU General Public License v2.0.
 * If a copy of the GPL was not distributed with this file, you can obtain one at:
 * http://www.gnu.org/licenses/gpl-2.0.html
 *
 */

#ifndef PHYSIKA_NUMERICS_LINEAR_SYSTEM_SOLVERS_LINEAR_SYSTEM_SOLVER_H_
#define PHYSIKA_NUMERICS_LINEAR_SYSTEM_SOLVERS_LINEAR_SYSTEM_SOLVER_H_

namespace Physika{

template <typename Scalar> class LinearSystem;
template <typename Scalar> class GeneralizedVector;

template <typename Scalar>
class LinearSystemSolver
{
public:
    LinearSystemSolver();
    virtual ~LinearSystemSolver();
    //procedure to solve the linear system
    bool solve(const LinearSystem<Scalar> &system, const GeneralizedVector<Scalar> &b, GeneralizedVector<Scalar> &x);
    //enable/disable preconditioner
    inline void enablePreconditioner() {use_preconditioner_ = true;}
    inline void disablePreconditioner() {use_preconditioner_ = false;}
    //query preconditioner enabled state
    inline bool isPreconditionerEnabled() const {return use_preconditioner_;}
protected:
    //virtual methods for subclass to implement
    virtual bool solveWithoutPreconditioner(const LinearSystem<Scalar> &system, const GeneralizedVector<Scalar> &b, GeneralizedVector<Scalar> &x) = 0;
    virtual bool solveWithPreconditioner(const LinearSystem<Scalar> &system, const GeneralizedVector<Scalar> &b, GeneralizedVector<Scalar> &x) = 0;
protected:
    bool use_preconditioner_;
};

}  //end of namespace Physika

#endif ///PHYSIKA_NUMERICS_LINEAR_SYSTEM_SOLVERS_LINEAR_SYSTEM_SOLVER_H_
