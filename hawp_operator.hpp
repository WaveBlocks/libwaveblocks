#ifndef WAVEBLOCKS_HAGEDORN_ABSTRACT_OPERATOR_HPP
#define WAVEBLOCKS_HAGEDORN_ABSTRACT_OPERATOR_HPP

class AbstractHaWp;
class AbstractHaWpOperator;

class AbstractHaWp
{
    friend class AbstractHaWpOperator;
    
public:
    virtual ~AbstractHaWp() { }
    
protected:
    virtual std::unique_ptr<AbstractHaWp> accept(AbstractHaWpOperator & op) const = 0;
};

class ScalarHaWp
{
protected:
    virtual std::unique_ptr<AbstractHaWp> accept(AbstractHaWpOperator & op) const
    {
        return op.visitLinearCombinedHaWp(*this);
    }
};

class LinearCombinedHaWp
{
public:
    std::vector<double> coefficients;
    std::vector< std::shared_ptr<AbstractHaWp> > vectors;
    
protected:
    virtual std::unique_ptr<AbstractHaWp> accept(AbstractHaWpOperator & op) const
    {
        return op.visitLinearCombinedHaWp(*this);
    }
};

class AbstractHaWpOperator
{
    friend class AbstractHaWp;
    
public:
    virtual ~AbstractHaWpOperator() { }
    
    std::unique_ptr<AbstractHaWp> operator()(AbstractHaWp const& wp)
    {
        return wp.accept(*this);
    }
    
protected:
    virtual std::unique_ptr<AbstractHaWp> visitScalarHaWp(ScalarHaWp const& wp) = 0;
    virtual std::unique_ptr<AbstractHaWp> visitLinearCombinedHaWp(LinearCombinedHaWp const& wp) = 0;
};

class GradientOperator : AbstractHaWpOperator
{
protected:
    virtual std::unique_ptr<AbstractHaWp> visitScalarHaWp(ScalarHaWp const& wp);
    virtual std::unique_ptr<AbstractHaWp> visitLinearCombinedHaWp(LinearCombinedHaWp const& wp);
};

#endif