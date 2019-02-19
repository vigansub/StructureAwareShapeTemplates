#ifndef GRAMMAR_H
#define GRAMMAR_H

#include "Headers.h"
#include "Connector.h"
#include <QDomDocument>
#include <time.h>


class Box;
class Node;

class Grammar : QObject
{
    Q_OBJECT
    
    enum Rand_Mode 
    {
        SINGLE = 0,
        CONSISTENT, 
        NONE, 
        GROUP
    };

	enum Independent3D
	{
		I3D_NONE = 0, 
		I3D_X,
		I3D_Y,
		I3D_XY,
		I3D_Z,
		I3D_XZ,
		I3D_YZ,
		I3D_XYZ
	};

	enum Independent2D
	{
		I2D_NONE = 0,
		I2D_X,
		I2D_Y,
		I2D_XY
	};

    enum Param_Type 
    {
        TYPE_NORMAL_DISTR_3D = 0, 
        TYPE_NORMAL_DISTR_2D,
        TYPE_RANGE_3D, 
		TYPE_RANGE_2D,
        TYPE_POSITION_2D, 
        TYPE_LOCATION, 
        TYPE_CONSTRAIN_SIZE
    };

    class Param
    {
    public:
        Param(const QString &_name, const QString &_setName, Param_Type _type, const QString &_mode) 
            : name(_name), type(_type), setName(_setName)
        {
            if(_mode == "consistent") { mode = CONSISTENT; }
            if(_mode == "single")     { mode = SINGLE; }
            if(_mode == "none")       { mode = NONE; }
            if(_mode == "group")      { mode = GROUP; }

            QString typeStr = "";
            if(type == TYPE_NORMAL_DISTR_3D)
                typeStr = "Normal_Distr_3D";
            if(type == TYPE_NORMAL_DISTR_2D)
                typeStr = "Normal_Distr_2D";
            if(type == TYPE_RANGE_3D)
                typeStr = "Range_3D";
            if(type == TYPE_RANGE_2D)
                typeStr = "Range_2D";
            if(type == TYPE_POSITION_2D)
                typeStr = "Position_2D";
            if(type == TYPE_LOCATION)
                typeStr = "Location";
            if(type == TYPE_CONSTRAIN_SIZE)
                typeStr = "Constrain_Size";

            globalId = _setName + "_" + typeStr + "_" + _name;
        }

        Param(const QString &_name, const QString &_globalId, Param_Type _type, Rand_Mode _mode) 
            : name(_name), type(_type), globalId(_globalId), mode(_mode) 
        {
        }

        virtual ~Param() {}

        QString name, globalId, setName;
        Param_Type type;
        Rand_Mode mode;
    };

    class Param3DNormalDistr : public Param
    {
    public:
        vec3 mu, sigma, propability, value, ddSigma, ddMu, resultMu, resultSigma;

        Param3DNormalDistr(const QString &_name, const QString &_setName, const QString &_mode, const vec3 &_mu, const vec3 &_sigma) 
            : Param(_name, _setName, TYPE_NORMAL_DISTR_3D, _mode), mu(_mu), sigma(_sigma), propability(1.0f, 1.0f, 1.0f), value(0.0f, 0.0f, 0.0f), ddSigma(0.0f, 0.0f, 0.0f), ddMu(0.0f, 0.0f, 0.0f), resultSigma(0.0f, 0.0f, 0.0f), resultMu(0.0f, 0.0f, 0.0f)
        {
            computeValueAndPropability();
            computeDerivatives();
        }

        Param3DNormalDistr(const QString &_name, const QString &_setName, const QString &_mode) 
            : Param(_name, _setName, TYPE_NORMAL_DISTR_3D, _mode), mu(0.0f, 0.0f, 0.0f), sigma(0.0f, 0.0f, 0.0f), propability(0.0f, 0.0f, 0.0f), value(0.0f, 0.0f, 0.0f), ddSigma(0.0f, 0.0f, 0.0f), ddMu(0.0f, 0.0f, 0.0f), resultSigma(0.0f, 0.0f, 0.0f), resultMu(0.0f, 0.0f, 0.0f)
        {
        }

        Param3DNormalDistr(const Param3DNormalDistr *p) : Param(p->name, p->globalId, p->type, p->mode)
        {
            this->mu = p->mu;
            this->sigma = p->sigma;
            this->propability = p->propability;
            this->value = p->value;
            this->ddSigma = p->ddSigma;
            this->ddMu = p->ddMu;
            this->resultMu = p->resultMu;
            this->resultSigma = p->resultSigma;
        }

        void computeValueAndPropability()
        {
            value.x = normalDistribution(mu.x, sigma.x);
            value.y = normalDistribution(mu.y, sigma.y);
            value.z = normalDistribution(mu.z, sigma.z);

            //if(value.x < 0.0f)
            //    value.x = 0.0f;
            //if(value.y < 0.0f)
            //    value.y = 0.0f;
            //if(value.z < 0.0f)
            //    value.z = 0.0f;

            propability.x = normalPDF(value.x, mu.x, sigma.x);
            propability.y = normalPDF(value.y, mu.y, sigma.y);           
            propability.z = normalPDF(value.z, mu.z, sigma.z);

            
        }

        void computeDerivatives()
        {
            //ddMu.x = normalDistributionDerivationMu(value.x, mu.x, sigma.x);
            //ddMu.y = normalDistributionDerivationMu(value.y, mu.y, sigma.y);
            //ddMu.z = normalDistributionDerivationMu(value.z, mu.z, sigma.z);

            //ddSigma.x = normalDistributionDerivationSigma(value.x, mu.x, sigma.x);
            //ddSigma.y = normalDistributionDerivationSigma(value.y, mu.y, sigma.y);
            //ddSigma.z = normalDistributionDerivationSigma(value.z, mu.z, sigma.z);

            ddMu.x = logOfNormalDistributionDerivationMu(value.x, mu.x, sigma.x);
            ddMu.y = logOfNormalDistributionDerivationMu(value.y, mu.y, sigma.y);
            ddMu.z = logOfNormalDistributionDerivationMu(value.z, mu.z, sigma.z);

            ddSigma.x = logOfNormalDistributionDerivationSigma(value.x, mu.x, sigma.x);
            ddSigma.y = logOfNormalDistributionDerivationSigma(value.y, mu.y, sigma.y);
            ddSigma.z = logOfNormalDistributionDerivationSigma(value.z, mu.z, sigma.z);
        }
    };

    class Param2DNormalDistr : public Param
    {
    public:
        vec2 mu, sigma, propability, value, ddSigma, ddMu, resultSigma, resultMu;

        Param2DNormalDistr(const QString &_name, const QString &_setName, const QString &_mode, const vec2 &_mu, const vec2 &_sigma) 
            : Param(_name, _setName, TYPE_NORMAL_DISTR_2D, _mode), mu(_mu), sigma(_sigma), propability(1.0f, 1.0f), value(0.0f, 0.0f), ddSigma(0.0f, 0.0f), ddMu(0.0f, 0.0f), resultSigma(0.0f, 0.0f), resultMu(0.0f, 0.0f)
        {
            computeValueAndPropability();
            computeDerivatives();
        }

        Param2DNormalDistr(const QString &_name, const QString &_setName, const QString &_mode) 
            : Param(_name, _setName, TYPE_NORMAL_DISTR_3D, _mode), mu(0.0f, 0.0f), sigma(0.0f, 0.0f), propability(0.0f, 0.0f), value(0.0f, 0.0f), ddSigma(0.0f, 0.0f), ddMu(0.0f, 0.0f), resultSigma(0.0f, 0.0f), resultMu(0.0f, 0.0f)
        {
        }

        Param2DNormalDistr(const Param2DNormalDistr *p) : Param(p->name, p->globalId, p->type, p->mode)
        {
            this->mu = p->mu;
            this->sigma = p->sigma;
            this->propability = p->propability;
            this->value = p->value;
            this->ddSigma = p->ddSigma;
            this->ddMu = p->ddMu;
            this->resultMu = p->resultMu;
            this->resultSigma = p->resultSigma;
        }


        void computeValueAndPropability()
        {
            value.x = normalDistribution(mu.x, sigma.x);            
            value.y = normalDistribution(mu.y, sigma.y);

            propability.x = normalPDF(value.x, mu.x, sigma.x);
            propability.y = normalPDF(value.y, mu.y, sigma.y);
        }

        void computeDerivatives()
        {
            //ddMu.x = normalDistributionDerivationMu(value.x, mu.x, sigma.x);
            //ddMu.y = normalDistributionDerivationMu(value.y, mu.y, sigma.y);

            //ddSigma.x = normalDistributionDerivationSigma(value.x, mu.x, sigma.x);
            //ddSigma.y = normalDistributionDerivationSigma(value.y, mu.y, sigma.y);

            ddMu.x = logOfNormalDistributionDerivationMu(value.x, mu.x, sigma.x);
            ddMu.y = logOfNormalDistributionDerivationMu(value.y, mu.y, sigma.y);

            ddSigma.x = logOfNormalDistributionDerivationSigma(value.x, mu.x, sigma.x);
            ddSigma.y = logOfNormalDistributionDerivationSigma(value.y, mu.y, sigma.y);
        }
    };

	class Param3DRange : public Param
	{
	public:
		vec3 m_mi, m_ma, m_value, propability;
		Independent3D independent;

		Param3DRange(const QString &_name, const QString &_setName, const QString &_mode, const vec3 &_mi, const vec3 &_ma)
			: Param(_name, _setName, TYPE_RANGE_3D, _mode), m_mi(_mi), m_ma(_ma), m_value(0.0f, 0.0f, 0.0f), propability(1.0f, 1.0f, 1.0f), independent(I3D_XYZ)
		{
			if (mode == CONSISTENT || mode == SINGLE || mode == NONE)
			{
				m_value.x = rand(m_mi.x, m_ma.x);
				m_value.y = rand(m_mi.y, m_ma.y);
				m_value.z = rand(m_mi.z, m_ma.z);
			}
		}

		Param3DRange(const QString &_name, const QString &_setName, const QString &_mode, const vec3 &_value)
			: Param(_name, _setName, TYPE_RANGE_3D, _mode), m_mi(), m_ma(), m_value(_value), independent(I3D_XYZ)
		{
		}

		void setIndependentVariables(Independent3D indep) { independent = indep; }
		Independent3D indpendentVariables() { return independent; }

		void setValue(const vec3 &v) { m_mi = v; m_ma = v; m_value = v; }
		void setValueIndep(const vec3 &v) 
		{
			//if (I3D_NONE) m_value = vec3(m_value.x, m_value.y, m_value.z);
			if (I3D_X) m_value = vec3(v.x, m_value.y, m_value.z);
			if (I3D_Y) m_value = vec3(m_value.x, v.y, m_value.z);
			if (I3D_XY) m_value = vec3(v.x, v.y, m_value.z);
			if (I3D_Z) m_value = vec3(m_value.x, m_value.y, v.z);
			if (I3D_XZ) m_value = vec3(v.x, m_value.y, v.z);
			if (I3D_YZ) m_value = vec3(m_value.x, v.y, v.z);
			if (I3D_XYZ) m_value = vec3(v.x, v.y, v.z);
		}

		vec3 value()
		{
			vec3 val = m_value;

			if (mode == SINGLE)
			{
				val.x = rand(m_mi.x, m_ma.x);
				val.y = rand(m_mi.y, m_ma.y);
				val.z = rand(m_mi.z, m_ma.z);
			}

			return val;
		}
    };

    class Param2DRange : public Param
    {
    public: 
        vec2 m_mi, m_ma, m_value, propability;
		Independent2D independent;

        Param2DRange(const QString &_name, const QString &_setName, const QString &_mode, const vec2 &_mi, const vec2 &_ma)
			: Param(_name, _setName, TYPE_RANGE_2D, _mode), m_mi(_mi), m_ma(_ma), m_value(0.0f, 0.0f), propability(1.0f, 1.0f), independent(I2D_XY)
        {
			if (mode == CONSISTENT || mode == SINGLE || mode == NONE)
			{
				m_value.x = rand(m_mi.x, m_ma.x);
				m_value.y = rand(m_mi.y, m_ma.y);
			}
        }

        Param2DRange(const QString &_name, const QString &_setName, const QString &_mode, const vec2 &_value)
			: Param(_name, _setName, TYPE_RANGE_3D, _mode), m_mi(), m_ma(), m_value(_value), independent(I2D_XY)
        {
        }

		void setIndependentVariables(Independent2D indep) { independent = indep; }
		Independent2D indpendentVariables() { return independent; }

		void setValue(const vec2 &v) { m_mi = v; m_ma = v; m_value = v; }
		void setValueIndep(const vec2 &v)
		{
			//if (I2D_NONE) m_value = vec2();
			if (I2D_X) m_value = vec2(v.x, m_value.y);
			if (I2D_Y) m_value = vec2(m_value.x, v.y);
			if (I2D_XY) m_value = vec2(v.x, v.y);
		}

		vec2 value()
		{
			vec2 val = m_value;

			if (mode == SINGLE)
			{
				val.x = rand(m_mi.x, m_ma.x);
				val.y = rand(m_mi.y, m_ma.y);
			}

			return val;
		}
    };

    class Param2DPosition: public Param
    {
    public:  
		vec2 m_value, m_offset;
		Independent2D independent;
        QString m_strX, m_strY;

        Param2DPosition(const QString &_name, const QString &_setName, const QString &_mode, const QString &_strX, const QString &_strY)
			: Param(_name, _setName, TYPE_POSITION_2D, _mode), m_value(0.0f, 0.0f), m_offset(0.0f, 0.0f), independent(I2D_NONE), m_strX(_strX), m_strY(_strY)
        {
			if (mode == NONE)
			{
				if (_strX == "center")
					m_value.x = 0.5f;
				if (_strX == "left")
					m_value.x = 0.0f;
				if (_strX == "right")
					m_value.x = 1.0f;
				if (_strX == "none")
					m_value.x = rand(0.0f, 1.0f);

				if (_strY == "center")
					m_value.y = 0.5f;
				if (_strY == "front")
					m_value.y = 0.0f;
				if (_strY == "back")
					m_value.y = 1.0f;
				if (_strY == "none") 
					m_value.y = rand(0.0f, 1.0f);
			}
            else if(mode == GROUP)
            {
				if (_strX == "center")
					m_value.x = 0.5f;
				if (_strX == "left")
					m_value.x = 0.0f + m_offset.x;
				if (_strX == "right")
					m_value.x = 1.0f - m_offset.x;

				if (_strY == "center")
					m_value.y = 0.5f;
				if (_strY == "front")
					m_value.y = 0.0f + m_offset.y;
				if (_strY == "back")
					m_value.y = 1.0f - m_offset.y;
            }
            else if (mode == CONSISTENT || mode == SINGLE)
            {
				m_value.x = rand(0.0f, 1.0f);
				m_value.y = rand(0.0f, 1.0f);

                //Something must change here to be able to enable independent variables
            }      
        }

		void setIndependentVariables(Independent2D indep) { independent = indep; }
		Independent2D indpendentVariables() { return independent; }

        void setOffset(const vec2 &v) { m_offset = v; }
		void setValue(const vec2 &v) { m_value = v; }
		void setValueIndep(const vec2 &v)
		{
			//if (I2D_NONE) m_value = vec2();
			if (I2D_X) m_value = vec2(v.x, m_value.y);
			if (I2D_Y) m_value = vec2(m_value.x, v.y);
			if (I2D_XY) m_value = vec2(v.x, v.y);
		}

        vec2 offset()
        {
            vec2 o(0.0f, 0.0f);

            if (mode == CONSISTENT || mode == SINGLE || mode == GROUP)
            {
			    if (m_strX == "left")
				    o.x = m_offset.x;
			    if (m_strX == "right")
				    o.x = -m_offset.x;

			    if (m_strY == "front")
				    o.y = m_offset.y;
			    if (m_strY == "back")
				    o.y = -m_offset.y;
            }
            return o;
        }

		vec2 value() 
		{
			vec2 val = m_value;
			if (mode == SINGLE)			
			{
				val.x = rand(0.0f, 1.0f);
				val.y = rand(0.0f, 1.0f);
				
			}

			return val;
		}
    };

    class ParamLocation : public Param
    {
    public:   
        Connector::Location location;
        float propability;

        ParamLocation(const QString &_name, const QString &_setName, const QString &_mode, const QString &type)
            : Param(_name, _setName, TYPE_LOCATION, _mode), propability(1.0f)
        {
            if(type == "bottom")
                location = Connector::BOTTOM;
            if(type == "top")
                location = Connector::TOP;
            if(type == "left")
                location = Connector::LEFT;
            if(type == "right")
                location = Connector::RIGHT;
            if(type == "front")
                location = Connector::FRONT;
            if(type == "back")
                location = Connector::BACK;     

            if(mode == CONSISTENT)
            {
                Connector::Location loc[] = {Connector::TOP, Connector::BOTTOM, Connector::RIGHT, Connector::LEFT, Connector::FRONT, Connector::BACK};        
                location = loc[rand() % 6];
                propability = 1.0f/6.0f;
            }
        }

        float propSumProd(bool sum = true)
        {
            return propability;
        }
    };

    class ParamConstrainSize: public Param
    {
    public:   
        int value;
        float propability;

        ParamConstrainSize(const QString &_name, const QString &_setName, const QString &_mode, int _value)
            : Param(_name, _setName, TYPE_CONSTRAIN_SIZE, _mode), value(_value), propability(1.0f)
        {
        }

    };

public:
    Grammar();
    Grammar(QString filename);
    ~Grammar();

    void init();
    void render(const Transform &trans, int shaderSelector = 0);
    void renderDepth(const Transform &trans, int shaderSelector = 0);
    void cleanUp();
    void softCleanUp();
    void update();

    void loadModelFromXML(QString fileName);
    void loadGrammarFromXML(QString fileName);
    void readParameter(const QString &setName, const QString &nodeName, const QDomElement &elem, vector<Param *> &params);
    void readOffsetSizeTransString(const QString &setName, const QDomElement &elem, vector<Param*> &params);

    void parse(QString &grammar);
    void derive();
#ifndef WIN32
    void parseTerminal(const QString &name, const QString &rule, vector<Param*> &overrideParams, int idx=-1);
    void parseTerminal(const QString &name, const QString &rule);
#else
    void parseTerminal(const QString &name, const QString &rule, vector<Param*> &overrideParams = vector<Param*>(), int idx=-1);
#endif
    void parseGroup(const QString &name, const QString &rule);
    QString selectRule(const QStringList &rules, const vector<float> &props);

	void adjustParameters(double trial[], int nr);
    void adjustParameters(const double* trial, int nr);
    vector<double> setParameterVals();
	void printParameters();
    
	void setIndependentVariablesNodes(vector<Param*> params, vector<Param*> paramsConnector, Node *n, Connector *c);

    void setupNode(Node *n, vector<Param*> &params);
    void setupConnector(Connector *c, vector<Param*> &params, vector<Param*> &overrideParams, int idx);
   
    vector<vec3> modelVertices();

	int determineNrIndependentVariables();
	QMultiMap<QString, Node*> nodes();  

public slots:
    void autoUpdate();
	
public:
	int m_nrIndependentVariables;

private:
    QMultiMap<QString, Node*> m_nodes;
    vector<Connector*> m_connectors; 

    Box *m_box;

    QTimer *m_timer;
    QDateTime m_oldTime;
    QString m_fileName;

    QMap<QString, pair<QStringList, vector<float>>> m_rules;
    QMap<QString, int> m_colorMap;

    //these maps hold the same data
    QMap<QString, vector<Param*>> m_parameterSets;
    QMap<QString, Param*> m_parameters;    
};

#endif

