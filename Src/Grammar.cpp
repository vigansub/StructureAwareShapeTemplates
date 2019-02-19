#include "Grammar.h"
#include "Box.h"
#include "Node.h"
#include "Connector.h"
#include "Headers.h"

#include <QDomDocument>
#include <QXmlStreamWriter>

#ifndef WIN32
Grammar::Grammar()
: m_nrIndependentVariables(0),
  m_fileName("./../_main/Data/Grammars/simple_exp2.gram")
#else
Grammar::Grammar()
: m_nrIndependentVariables(0), m_fileName("Data/Grammars/chair/chair_gram3.gram")
#endif
{
    m_timer = new QTimer();
    connect(m_timer, SIGNAL(timeout()), this, SLOT(autoUpdate()));  
    m_timer->start(1000);    

    m_box = new Box();

    autoUpdate();

    //loadGrammarFromXML(m_fileName);
    //loadModelFromXML(m_fileName);
    //init();

}

#ifndef WIN32
Grammar::Grammar(QString filename)
    : m_nrIndependentVariables(0)
    #else
Grammar::Grammar(QString filename)
    : m_nrIndependentVariables(0)
    #endif
{
    m_fileName = filename;
    m_timer = new QTimer();
    connect(m_timer, SIGNAL(timeout()), this, SLOT(autoUpdate()));
    m_timer->start(1000);

    m_box = new Box();

    autoUpdate();

    //loadGrammarFromXML(m_fileName);
    //loadModelFromXML(m_fileName);
    //init();

}

Grammar::~Grammar()
{
}

void Grammar::init()
{    
    Node *n1 = new Node(m_box, "1");
    Node *n2 = new Node(m_box, "2");
    Node *n3 = new Node(m_box, "3");
    Node *n4 = new Node(m_box, "4");
    Node *n5 = new Node(m_box, "5");
    
    n4->setSize(vec3(4, 4, 4));

    Connector *c1 = new Connector(m_box, n1, n2);
    c1->setLocation(Connector::RIGHT);
    c1->setSize(0.1, 0.1);
    c1->setOffset(0.5, 1.0);

    n1->addChild(c1);
    n2->addParent(c1);

    Connector *c2 = new Connector(m_box, n2, n3);
    c2->setLocation(Connector::FRONT);
    c2->setSize(0.5, 0.5);
    c2->setOffset(1.0, 0);

    n2->addChild(c2);
    n3->addParent(c2);

    Connector *c3 = new Connector(m_box, n1, n4);
    c3->setLocation(Connector::TOP);
    c3->setSize(0.5, 0.5);
    c3->setOffset(0.0, 1);
    
    n1->addChild(c3);
    n4->addParent(c3);

    Connector *c4 = new Connector(m_box, n4, n5);
    c4->setLocation(Connector::LEFT);
    c4->setSize(0.5, 0.5);
    c4->setOffset(0.0, 1);
    
    n4->addChild(c4);
    n5->addParent(c4);


    m_nodes.insert(n1->name(), n1);
    m_nodes.insert(n2->name(), n2);
    m_nodes.insert(n3->name(), n3);
    m_nodes.insert(n4->name(), n4);
    m_nodes.insert(n5->name(), n5);

    m_connectors.push_back(c1);
    m_connectors.push_back(c2);
    m_connectors.push_back(c3);
    m_connectors.push_back(c4);
}

void Grammar::render(const Transform &trans, int shaderSelector)
{
    if(params::inst()->renderConnectors)
    {
        for(int i=0; i<m_connectors.size(); ++i)
        {
            m_connectors[i]->render(trans, shaderSelector);
        }
    }

    if(params::inst()->renderBoxes)
    {
        for(auto iter=m_nodes.begin(); iter!=m_nodes.end(); ++iter)
        {
         iter.value()->render(trans, shaderSelector);
        }
    }
}

void Grammar::renderDepth(const Transform &trans, int shaderSelector)
{
    for (auto iter = m_nodes.begin(); iter != m_nodes.end(); ++iter)
    {
        iter.value()->renderDepth(trans, shaderSelector);
    }
}

void Grammar::loadModelFromXML(QString fileName)
{
    QFile file(fileName);

	if (!file.open(QIODevice::ReadOnly))
		return;

	QDomDocument doc(fileName);
	if (!doc.setContent(&file))        
		return;

	file.close();

	QDomNodeList listShapes = doc.elementsByTagName("Shape");
	QDomNode     nodeShapes = listShapes.at(0);
	QDomNodeList nodesShapes = nodeShapes.childNodes();

	for(int i=0; i<nodesShapes.size(); ++i)
    {
		QDomNode node = nodesShapes.at(i);

        if(node.nodeName().compare("Part") == 0)
        {
		    QDomNodeList nodeTags = nodesShapes.at(i).childNodes();

            QDomElement element = node.toElement();
		    QString name = element.attribute("name");
            QString conSize = element.attribute("constrainSize");            

            bool constrainSize = false;
            if(conSize.compare("true") == 0)
                constrainSize = true;           

            Node *part = new Node(m_box, name);

		    for(int j=0; j<nodeTags.size(); ++j)
            {
			    QDomNode n = nodeTags.at(j);
			    QDomElement e = n.toElement();

  			    if(n.nodeName().compare("Size") == 0)
			    {
				    QString strX = e.attribute("x");                
				    QString strY = e.attribute("y");
				    QString strZ = e.attribute("z");

                    vec3 size = vec3(strX.toFloat(), strY.toFloat(), strZ.toFloat());
                    part->setSize(size);                    
			    }			

  			    if(n.nodeName().compare("Translation") == 0)
			    {
				    QString strX = e.attribute("x");                
				    QString strY = e.attribute("y");
				    QString strZ = e.attribute("z");

                    vec3 trans = vec3(strX.toFloat(), strY.toFloat(), strZ.toFloat());
                    part->setTranslation(trans); 
			    }			

		    }

            m_nodes.insert(part->name(), part);
        }


        if(node.nodeName().compare("Connector") == 0)
        {
		    QDomNodeList nodeTags = nodesShapes.at(i).childNodes();

            QDomElement element = node.toElement();
		    
            QString strFrom = element.attribute("from");
            QString strTo = element.attribute("to");

            Node *from = m_nodes.find(strFrom).value();
            Node *to = m_nodes.find(strTo).value();

            QString typeStr = element.attribute("type");
            Connector::Location location;

            if(typeStr.compare("top") == 0)
                location = Connector::TOP;
            if(typeStr.compare("bottom") == 0)
                location = Connector::BOTTOM;
            if(typeStr.compare("left") == 0)
                location = Connector::LEFT;
            if(typeStr.compare("right") == 0)
                location = Connector::RIGHT;
            if(typeStr.compare("front") == 0)
                location = Connector::FRONT;
            if(typeStr.compare("back") == 0)
                location = Connector::BACK;

            Connector *c = new Connector(m_box, from, to);
            c->setLocation(location);

		    for(int j=0; j<nodeTags.size(); ++j)
            {
			    QDomNode n = nodeTags.at(j);
			    QDomElement e = n.toElement();

  			    if(n.nodeName().compare("Size") == 0)
			    {
				    QString strX = e.attribute("x");                
				    QString strY = e.attribute("y");

                    c->setSize(strX.toFloat(), strY.toFloat());
			    }			

  			    if(n.nodeName().compare("Offset") == 0)
			    {
				    QString strX = e.attribute("x");                
				    QString strY = e.attribute("y");

                    c->setOffset(strX.toFloat(), strY.toFloat());
			    }			

		    }

            from->addChild(c);
            to->addParent(c);

            m_connectors.push_back(c);
        }
    }
}

void Grammar::loadGrammarFromXML(QString fileName)
{
    QFile file(fileName);

	if (!file.open(QIODevice::ReadOnly))
		return;

	QDomDocument doc(fileName);
	if (!doc.setContent(&file))        
		return;

	file.close();

	QDomNodeList listShapes = doc.elementsByTagName("Shape");
	QDomNode     nodeShapes = listShapes.at(0);
	QDomNodeList nodesShapes = nodeShapes.childNodes();

	for(int i=0; i<nodesShapes.size(); ++i)
    {
		QDomNode node = nodesShapes.at(i);

        if(node.nodeName().compare("Grammar") == 0)
        {
		    QDomNodeList nodeTags = nodesShapes.at(i).childNodes();
            QDomElement element = node.toElement();		    
            QString text = element.text();

            parse(text);
        }

        if(node.nodeName().compare("Parameters") == 0)
        {
            QDomNodeList nodeTags = nodesShapes.at(i).childNodes();
            QDomElement element = node.toElement();
		    
            int colorCounter = 0;

		    for(int j=0; j<nodeTags.size(); ++j)
            {
                QDomNode n = nodeTags.at(j);
			    QDomElement e = n.toElement();
                QDomNodeList childTags = n.childNodes();

                QString setName = n.nodeName();
                QString type = e.attribute("type");
                
                m_colorMap.insert(setName, colorCounter);
                colorCounter++;

                vector<Param*> params;

                for(int k=0; k<childTags.size(); ++k)
                {
			        QDomNode c = childTags.at(k);
			        QDomElement ec = c.toElement();

                    readParameter(setName, c.nodeName(), ec, params);
                }

                m_parameterSets.insert(setName, params);
		    }
        }
    }
}

void Grammar::parse(QString &grammar)
{
    QString g = grammar.replace(" ", "");
    QStringList rulesStr = g.split("\r");

    for(int i=0; i<rulesStr.size(); ++i)
    {
        QString cleanRule = rulesStr[i].replace("\n", "");
        QStringList leftRight = cleanRule.split("-");
        
        if(leftRight.size() > 1)
        {
            pair<QStringList, vector<float>> rulesProp;

            if(leftRight[1].contains("{") && leftRight[1].contains("}"))
            {
                QString tmp = leftRight[1];
                tmp.chop(1);

                QStringList list = tmp.split("{");

                QString left = list[0];
                QString right = list[1];

                QStringList rules = left.split("|");
                QStringList tmpProps = right.split(",");

                vector<float> props;
                for(int j=0; j<tmpProps.size(); ++j)
                {
                    props.push_back(tmpProps[j].toFloat());
                }

                rulesProp.first = rules;
                rulesProp.second = props;
            }
            else
            {
                QStringList rules = leftRight[1].split("|");
                vector<float> props;

                rulesProp.first = rules;
                rulesProp.second = props;
            }          
			
            m_rules.insert(leftRight[0], rulesProp);
        }
    }
}

void Grammar::derive()
{
    //srand(time(0));

    QString d = "R";

    int idx = 0;

    while(idx < d.length())
    {
        QString s(d[idx]);
        auto iter = m_rules.find(s);
        
        //random
        pair<QStringList, vector<float>> &entry = iter.value();
        QString r = selectRule(entry.first, entry.second);
     
        if(r.contains("n") || r.contains("c"))
        {
            parseTerminal(s, r);
            d = d.remove(idx, 1);
            idx = 0;
        }
        else if(r.contains("[") && r.contains("]"))
        {
            parseGroup(s, r);
            d = d.remove(idx, 1);
            idx = 0;
        }
        else if(r.contains("0"))
        {
            d = d.remove(idx, 1);
            idx = 0;
        }
        else
        {
            d.replace(s, r);
            idx = 0;
        }                
    }
}

QString Grammar::selectRule(const QStringList &rules, const vector<float> &props)
{
    QString s;

    if(props.size() > 0)
    {
        //random
        float rnd = rand(0.0, 1.0);        

        float lb = 0.0f;
        float ub = 0.0f;               

        for(int i=0; i<rules.size(); ++i)
        {
            float c = props[i];
            ub = lb + c;

            if(rnd >= lb && rnd < ub)
            {                
                s = rules[i];
                
                //not implemented right now
                //m_propDens.push_back(PropInfo(s, props[i]));

                break;
            }

            lb = ub;
        }   
    }
    else
    {
        s = rules[rand()%rules.size()];
        //not implemented right now
        //m_propDens.push_back(PropInfo(s, 1.0f/rules.size()));
    }

    return s;
}

void Grammar::parseGroup(const QString &name, const QString &rule)
{
    QString setName = QString(rule[rule.length()-3]) + QString(rule[rule.length()-2]);

    QString tmp = rule;
    tmp.chop(4);
    tmp = tmp.remove("[");
    tmp = tmp.remove("]");    

    vector<Param*> &pi = m_parameterSets.find(setName).value();
    vector<Param*> piFinal(pi.begin()+1, pi.end());

    //for(auto iter = m_parameters.begin(); iter!=m_parameters.end(); ++iter)
    //    qDebug() << iter.value()->globalId;

    //qDebug() << setName << rule << name << tmp;

    QString paramStr = setName + "_Range_2D_Offset";
    auto paramEntry = m_parameters.find(paramStr);
    
    vec2 offset(0.0f, 0.0f);
    if(paramEntry != m_parameters.end())
    {
        Param2DRange *p = dynamic_cast<Param2DRange*>(paramEntry.value());
        offset = p->value();
    }

    for(int i=0; i<piFinal.size(); ++i)
    {
        if(piFinal[i]->name.contains("Position"))
        {
            Param2DPosition *t = dynamic_cast<Param2DPosition *>(piFinal[i]);
            t->setOffset(offset);
        }
    }

    for(int i=0; i<tmp.length(); ++i)
    {
        QString s(tmp[i]);
        auto iter = m_rules.find(s);        

        if(iter != m_rules.end())
        {
            int rnd = 0;//rand()%iter.value().size();            
            QString r = iter.value().first[rnd];

            parseTerminal(s, r, piFinal, i);
        }
    }
}

void Grammar::parseTerminal(const QString &name, const QString &rule, vector<Param*> &overrideParams, int idx)
{
    //Node and Connector
    if(rule.length() % 12 == 0)
    {
        for(int i=0; i<rule.length(); i+=12)
        {
            //Node
            Node *n = new Node(m_box, name);        
            vector<Param*> &paramsNode = m_parameterSets[QString(rule[i+2]) + QString(rule[i+3])];
            int colorMode = m_colorMap[QString(rule[i+2]) + QString(rule[i+3])];
            n->setColorMode(colorMode);

            setupNode(n, paramsNode);

            
            //Connector
            Node *parent = m_nodes.find(QString(rule[i+7])).value();
            Connector * c = new Connector(m_box, parent, n);

            vector<Param*> &paramsConnector = m_parameterSets[QString(rule[i+9]) + QString(rule[i+10])];
            setupConnector(c, paramsConnector, overrideParams, idx);
                
            if(parent->addChild(c))
            {
                n->addParent(c);
            }

			setIndependentVariablesNodes(paramsNode, paramsConnector, n, c);

            m_nodes.insertMulti(name, n);
            m_connectors.push_back(c);
        }
    }

    //Node only
    if(rule.length() == 5)
    {
        //Node
        Node *n = new Node(m_box, name);        
        vector<Param*> &paramInfo = m_parameterSets[QString(rule[2]) + QString(rule[3])];
        setupNode(n, paramInfo);

        m_nodes.insertMulti(name, n);
    }
}

#ifndef WIN32
void Grammar::parseTerminal(const QString &name, const QString &rule)
{
    vector<Param*> temp;
    parseTerminal(name, rule, temp);
}
#endif

void Grammar::setIndependentVariablesNodes(vector<Param*> params, vector<Param*> paramsConnector, Node *n, Connector *c)
{
    int cMode = 0;

    for(int i=0; i<params.size(); ++i)
    {
		Param *p = params[i];
		if (p->type == TYPE_CONSTRAIN_SIZE)
		{
            ParamConstrainSize *s = dynamic_cast<ParamConstrainSize*> (p);
            cMode = s->value;
        }
    }

    Param2DRange *s_connect;

    //qDebug() << paramsConnector.size();

    for (int i = 0; i < paramsConnector.size(); ++i)
    {
        Param *p = paramsConnector[i];
        //qDebug() << p->globalId;
        if (p->type == TYPE_RANGE_2D)
        {
            s_connect = dynamic_cast<Param2DRange *> (p);
            //qDebug() << s_connect->independent;
        }
    }

	if (c->location() == Connector::TOP)
	{
		for (int j = 0; j < params.size(); ++j)
		{
			Param *p = params[j];

            if (p->type == TYPE_RANGE_3D)
            {
                Param3DRange *s = dynamic_cast<Param3DRange *> (p);

                if (cMode == 0)
                    s->setIndependentVariables(I3D_Y);
                else if (cMode == 1)
                {
                    s->setIndependentVariables(I3D_YZ);
                    s_connect->setIndependentVariables(I2D_X);
                }
                else if (cMode == 2)
                {
                    s->setIndependentVariables(I3D_XY);
                    s_connect->setIndependentVariables(I2D_Y);
                }
                else if (cMode == 3)
                {
                    s->setIndependentVariables(I3D_XYZ);
                    s_connect->setIndependentVariables(I2D_NONE);
                }
            }
		}
	}

	if (c->location() == Connector::BOTTOM)
	{
		for (int j = 0; j < params.size(); ++j)
		{
			Param *p = params[j];

            if (p->type == TYPE_RANGE_3D)
            {
                Param3DRange *s = dynamic_cast<Param3DRange *> (p);

                if (cMode == 0)
                    s->setIndependentVariables(I3D_Y);
                else if (cMode == 1)
                {
                    s->setIndependentVariables(I3D_YZ);
                    s_connect->setIndependentVariables(I2D_X);
                }
                else if (cMode == 2)
                {
                    s->setIndependentVariables(I3D_XY);
                    s_connect->setIndependentVariables(I2D_Y);
                }
                else if (cMode == 3)
                {
                    s->setIndependentVariables(I3D_XYZ);
                    s_connect->setIndependentVariables(I2D_NONE);
                }
            }
		}
	}

	if (c->location() == Connector::FRONT)
	{
		for (int j = 0; j < params.size(); ++j)
		{
			Param *p = params[j];

            if (p->type == TYPE_RANGE_3D)
            {
                Param3DRange *s = dynamic_cast<Param3DRange *> (p);

                if (cMode == 0)
                    s->setIndependentVariables(I3D_Z);
                else if (cMode == 1)
                {
                    s->setIndependentVariables(I3D_YZ);
                    s_connect->setIndependentVariables(I2D_X);
                }
                else if (cMode == 2)
                {
                    s->setIndependentVariables(I3D_XZ);
                    s_connect->setIndependentVariables(I2D_Y);
                }
                else if (cMode == 3)
                {
                    s->setIndependentVariables(I3D_XYZ);
                    s_connect->setIndependentVariables(I2D_NONE);
                }
            }
		}
	}

	if (c->location() == Connector::BACK)
	{
		for (int j = 0; j < params.size(); ++j)
		{
			Param *p = params[j];

            if (p->type == TYPE_RANGE_3D)
            {
                Param3DRange *s = dynamic_cast<Param3DRange *> (p);

                if (cMode == 0)
                    s->setIndependentVariables(I3D_Z);
                else if (cMode == 1)
                {
                    s->setIndependentVariables(I3D_YZ);
                    s_connect->setIndependentVariables(I2D_X);
                }
                else if (cMode == 2)
                {
                    s->setIndependentVariables(I3D_XZ);
                    s_connect->setIndependentVariables(I2D_Y);
                }
                else if (cMode == 3)
                {
                    s->setIndependentVariables(I3D_XYZ);
                    s_connect->setIndependentVariables(I2D_NONE);
                }
            }
		}
	}

	if (c->location() == Connector::RIGHT)
	{
		for (int j = 0; j < params.size(); ++j)
		{
			Param *p = params[j];

            if (p->type == TYPE_RANGE_3D)
            {
                Param3DRange *s = dynamic_cast<Param3DRange *> (p);

                if (cMode == 0)
                    s->setIndependentVariables(I3D_X);
                else if (cMode == 1)
                {
                    s->setIndependentVariables(I3D_XY);
                    s_connect->setIndependentVariables(I2D_X);
                }
                else if (cMode == 2)
                {
                    s->setIndependentVariables(I3D_XZ);
                    s_connect->setIndependentVariables(I2D_Y);
                }
                else if (cMode == 3)
                {
                    s->setIndependentVariables(I3D_XYZ);
                    s_connect->setIndependentVariables(I2D_NONE);
                }
            }
		}
	}

	if (c->location() == Connector::LEFT)
	{
		for (int j = 0; j < params.size(); ++j)
		{
			Param *p = params[j];

            if (p->type == TYPE_RANGE_3D)
            {
                Param3DRange *s = dynamic_cast<Param3DRange *> (p);

                if (cMode == 0)
                    s->setIndependentVariables(I3D_X);
                else if (cMode == 1)
                {
                    s->setIndependentVariables(I3D_XY);
                    s_connect->setIndependentVariables(I2D_X);
                }
                else if (cMode == 2)
                {
                    s->setIndependentVariables(I3D_XZ);
                    s_connect->setIndependentVariables(I2D_Y);
                }
                else if (cMode == 3)
                {
                    s->setIndependentVariables(I3D_XYZ);
                    s_connect->setIndependentVariables(I2D_NONE);
                }
            }
		}
	}
}

void Grammar::setupNode(Node *n, vector<Param*> &params)
{
    for(int i=0; i<params.size(); ++i)
    {
        Param *p = params[i];        

        if(p->name == "Translation")
        {
            Param3DRange *t = dynamic_cast<Param3DRange *>(p);
            n->setTranslation(t->value());
        }

        if(p->name == "Constrain")
        {
            ParamConstrainSize *t = dynamic_cast<ParamConstrainSize *>(p);

            n->setConstrainSize(t->value);
        }

        if(p->name == "Size")
        {
            Param3DRange *t = dynamic_cast<Param3DRange *>(p);
            n->setSize(t->value());
        }
    }
}

void Grammar::setupConnector(Connector *c, vector<Param*> &params,  vector<Param*> &overrideParams, int idx)
{
    for(int j=0; j<params.size(); ++j)
    {
        Param *p = params[j];

        ////offset
        //if(p.name == "Offset")
        //{
        //    vec3 mi = p.mi; 
        //    vec3 ma = p.ma;
        //    if(idx >= 0)
        //    {
        //        mi = overrideParams[idx].mi;
        //        ma = overrideParams[idx].ma;
        //    }
        //    //c->setOffsetRnd(mi, ma);
        //    c->setOffset(0.5, 0.5);

        //}

        if(p->name == "Position")
        {
            Param2DPosition *t = dynamic_cast<Param2DPosition *>(p);
            vec2 pos = t->value(); 
            vec2 offset = vec2(0.0f, 0.0f);

            if(idx >= 0)
            {
                Param2DPosition *s = dynamic_cast<Param2DPosition *> (overrideParams[idx]);
                pos = s->value();
                offset = s->offset();
            }

            c->adjustOffset(pos.x + offset.x, pos.y + offset.y);
        }

        if(p->name == "Size")
        {
            Param2DRange *t = dynamic_cast<Param2DRange *>(p);
			vec2 val = t->value();

            c->setSize(val.x, val.y);
            
        }

        if(p->name == "Location")
        {
            ParamLocation *t = dynamic_cast<ParamLocation *>(p);
            c->setLocation(t->location);
        }
    }
}

void Grammar::cleanUp()
{
    softCleanUp();
    
    for(auto iter = m_parameters.begin(); iter!=m_parameters.end(); ++iter)
    {
        delete iter.value();
    }

    m_parameters.clear();
    m_rules.clear();    
    m_parameterSets.clear();
}

void Grammar::softCleanUp()
{
    for(auto iter = m_nodes.begin(); iter!=m_nodes.end(); ++iter)
    {
        Node *n = iter.value();
        delete n;
    }    

    for(int i=m_connectors.size()-1; i>=0; --i)
    {
        Connector *c = m_connectors[i];
        delete c;
    }

    m_nodes.clear();        
    m_connectors.clear();
}

void Grammar::update()
{
    cleanUp();    
    loadGrammarFromXML(m_fileName);
    derive();

    //Output Parameters
    //for(auto iter = m_parameters.begin(); iter!=m_parameters.end(); ++iter)
    //{
    //    qDebug() << iter.key();
    //}

    //for(auto iter = m_propDens.begin(); iter!=m_propDens.end(); ++iter)
    //{
    //    qDebug() << iter.key() << iter.value();
    //}
}

void Grammar::autoUpdate()
{
    QString fName(m_fileName);
    QFileInfo fileInfo(fName);
    QDateTime dt = fileInfo.lastModified();  

    if(dt != m_oldTime)
    {
        update();
        m_oldTime = dt;

		for (auto iter = m_parameters.begin(); iter != m_parameters.end(); ++iter)
		{
			//qDebug() << iter.key();
		}

        qDebug() << "IndependentVariables: " << determineNrIndependentVariables();
    }    
}

void Grammar::readParameter(const QString &setName, const QString &nodeName, const QDomElement &elem, vector<Param *> &params)
{
    QString name = elem.nodeName();
    QString rand_mode = elem.attribute("rand_mode");
    
    if(nodeName == "Translation")
        readOffsetSizeTransString(setName, elem, params);

  	if(nodeName == "Size")
        readOffsetSizeTransString(setName, elem, params);

  	if(nodeName.contains("Offset"))
        readOffsetSizeTransString(setName, elem, params);

  	if(nodeName.contains("Position"))
    {
        Param *p = new Param2DPosition(name, setName, rand_mode, elem.attribute("x"), elem.attribute("y"));
        params.push_back(p);
        m_parameters.insert(p->globalId, p);
    }

  	if(nodeName == "Location")
    {
        QString type = elem.attribute("type");    

        Param *p = new ParamLocation(name, setName, rand_mode, type);
        params.push_back(p);
        m_parameters.insert(p->globalId, p);
    }

  	if(nodeName == "Constrain")
    {
        Param *p = new ParamConstrainSize(name, setName, rand_mode, elem.attribute("value").toInt());
        params.push_back(p);
        m_parameters.insert(p->globalId, p);
    }
}

void Grammar::readOffsetSizeTransString(const QString &setName, const QDomElement &elem, vector<Param*> &params)
{
    if(elem.hasAttribute("xMu") && elem.hasAttribute("xSigma"))
    {
        QString name = elem.nodeName();
        QString rand_mode = elem.attribute("rand_mode");

        QString strXSigma = elem.attribute("xSigma");
        QString strXMu = elem.attribute("xMu");

        QString strYSigma = elem.attribute("ySigma");
        QString strYMu = elem.attribute("yMu");

        if(elem.hasAttribute("zMu") && elem.hasAttribute("zSigma"))
        {
            QString strZSigma = elem.attribute("zSigma");
            QString strZMu = elem.attribute("zMu");

            vec3 sigma = vec3(strXSigma.toFloat(), strYSigma.toFloat(), strZSigma.toFloat());
            vec3 mu = vec3(strXMu.toFloat(), strYMu.toFloat(), strZMu.toFloat());

            Param *p = new Param3DNormalDistr(name, setName, rand_mode, mu, sigma);
            params.push_back(p);
            m_parameters.insert(p->globalId, p);
        }
        else
        {
            vec2 sigma = vec2(strXSigma.toFloat(), strYSigma.toFloat());
            vec2 mu = vec2(strXMu.toFloat(), strYMu.toFloat());

            Param *p = new Param2DNormalDistr(name, setName, rand_mode, mu, sigma);
            params.push_back(p);
            m_parameters.insert(p->globalId, p);
        }
    }
    else
    {
        QString name = elem.nodeName();
        QString rand_mode = elem.attribute("rand_mode");

        QString strXFrom = elem.attribute("xfrom");
        QString strXTo = elem.attribute("xto");

        QString strYFrom = elem.attribute("yfrom");
        QString strYTo = elem.attribute("yto");

        if(elem.hasAttribute("zfrom") && elem.hasAttribute("zto"))
        {
            QString strZFrom = elem.attribute("zfrom");
            QString strZTo = elem.attribute("zto");

            vec3 from = vec3(strXFrom.toFloat(), strYFrom.toFloat(), strZFrom.toFloat());
            vec3 to = vec3(strXTo.toFloat(), strYTo.toFloat(), strZTo.toFloat());

            Param *p = new Param3DRange(name, setName, rand_mode, from, to);
            params.push_back(p);
            m_parameters.insert(p->globalId, p);
        }
        else
        {
            vec2 from = vec2(strXFrom.toFloat(), strYFrom.toFloat());
            vec2 to = vec2(strXTo.toFloat(), strYTo.toFloat());

            Param *p = new Param2DRange(name, setName, rand_mode, from, to);
            params.push_back(p);
            m_parameters.insert(p->globalId, p);
        }
    }
}

vector<vec3> Grammar::modelVertices()
{
    vector<vec3> modelVertices;
    for(auto iter=m_nodes.begin(); iter!=m_nodes.end(); ++iter)
    {
        Node *n = iter.value();
        
        mat4 model = n->tranform();
        vector<Vertex> verts = m_box->vertices();

        for(int i=0; i<verts.size(); ++i)
        {
            vec3 v = verts[i].position;
            vec4 t = model * vec4(v.x, v.y, v.z, 1.0f);
            v = vec3(t.x, t.y, t.z);

            modelVertices.push_back(v);
        }
    }

    return modelVertices;
}

void Grammar::adjustParameters(double trial[], int nr)
{        
    //softCleanUp();
    //derive();

	int idx = 0;
	for (auto iter = m_parameters.begin(); iter != m_parameters.end(); ++iter)
	{
		Param *p = iter.value();

		if (p->type == TYPE_RANGE_3D)
		{
			Param3DRange *t = dynamic_cast<Param3DRange *>(p);
			Independent3D indep = t->indpendentVariables();

			vec3 val = t->value();
			
			if (indep == I3D_X)
			{
				val.x = trial[idx];
				idx++;
			}
			if (indep == I3D_Y)
			{
				val.y = trial[idx];
				idx++;
			}
			if (indep == I3D_Z)
			{
				val.z = trial[idx];
				idx++;
			}
			if (indep == I3D_XY)
			{
				val.x = trial[idx];
				val.y = trial[idx + 1];
				idx += 2;
			}
			if (indep == I3D_YZ)
			{
				val.y = trial[idx];
				val.z = trial[idx + 1];
				idx += 2;
			}
			if (indep == I3D_XZ)
			{
				val.x = trial[idx];
				val.z = trial[idx + 1];
				idx += 2;
			}
			if (indep == I3D_XYZ)
			{
				val.x = trial[idx];
				val.y = trial[idx + 1];
				val.z = trial[idx + 2];
				idx += 3;
			}

			t->setValue(val);
		}

		if (p->type == TYPE_RANGE_2D)
		{
			Param2DRange *t = dynamic_cast<Param2DRange *>(p);
			Independent2D indep = t->indpendentVariables();

			vec2 val = t->value();

			if (indep == I2D_X)
			{
				val.x = trial[idx];
				idx++;
			}
			if (indep == I2D_Y)
			{
				val.y = trial[idx];
				idx++;
			}
			if (indep == I2D_XY)
			{
				val.x = trial[idx];
				val.y = trial[idx + 1];
				idx += 2;
			}

			t->setValue(val);
		}

		if (p->type == TYPE_POSITION_2D)
		{
			Param2DPosition *t = dynamic_cast<Param2DPosition *>(p);
			Independent2D indep = t->indpendentVariables();

			vec2 val = t->value();

			if (indep == I2D_X)
			{
				val.x = trial[idx];
				idx++;
			}
			if (indep == I2D_Y)
			{
				val.y = trial[idx];
				idx++;
			}
			if (indep == I2D_XY)
			{
				val.x = trial[idx];
				val.y = trial[idx + 1];
				idx += 2;
			}

			t->setValue(val);
		}
	}
}

void Grammar::adjustParameters(const double* trial, int nr)
{
    //softCleanUp();
    //derive();

    int idx = 0;
    for (auto iter = m_parameters.begin(); iter != m_parameters.end(); ++iter)
    {
        Param *p = iter.value();

        if (p->type == TYPE_RANGE_3D)
        {
            Param3DRange *t = dynamic_cast<Param3DRange *>(p);
            Independent3D indep = t->indpendentVariables();

            vec3 val = t->value();

            if (indep == I3D_X)
            {
                val.x = trial[idx];
                idx++;
            }
            if (indep == I3D_Y)
            {
                val.y = trial[idx];
                idx++;
            }
            if (indep == I3D_Z)
            {
                val.z = trial[idx];
                idx++;
            }
            if (indep == I3D_XY)
            {
                val.x = trial[idx];
                val.y = trial[idx + 1];
                idx += 2;
            }
            if (indep == I3D_YZ)
            {
                val.y = trial[idx];
                val.z = trial[idx + 1];
                idx += 2;
            }
            if (indep == I3D_XZ)
            {
                val.x = trial[idx];
                val.z = trial[idx + 1];
                idx += 2;
            }
            if (indep == I3D_XYZ)
            {
                val.x = trial[idx];
                val.y = trial[idx + 1];
                val.z = trial[idx + 2];
                idx += 3;
            }

            t->setValue(val);
        }

        if (p->type == TYPE_RANGE_2D)
        {
            Param2DRange *t = dynamic_cast<Param2DRange *>(p);
            Independent2D indep = t->indpendentVariables();

            vec2 val = t->value();

            if (indep == I2D_X)
            {
                val.x = trial[idx];
                idx++;
            }
            if (indep == I2D_Y)
            {
                val.y = trial[idx];
                idx++;
            }
            if (indep == I2D_XY)
            {
                val.x = trial[idx];
                val.y = trial[idx + 1];
                idx += 2;
            }

            t->setValue(val);
        }

        if (p->type == TYPE_POSITION_2D)
        {
            Param2DPosition *t = dynamic_cast<Param2DPosition *>(p);
            Independent2D indep = t->indpendentVariables();

            vec2 val = t->value();

            if (indep == I2D_X)
            {
                val.x = trial[idx];
                idx++;
            }
            if (indep == I2D_Y)
            {
                val.y = trial[idx];
                idx++;
            }
            if (indep == I2D_XY)
            {
                val.x = trial[idx];
                val.y = trial[idx + 1];
                idx += 2;
            }

            t->setValue(val);
        }
    }
}

vector<double> Grammar::setParameterVals()
{
    //softCleanUp();
    //derive();

    int idx = 0;
    vector<double> all_params;
    all_params.resize(determineNrIndependentVariables());
    for (auto iter = m_parameters.begin(); iter != m_parameters.end(); ++iter)
    {
        Param *p = iter.value();

        if (p->type == TYPE_RANGE_3D)
        {
            Param3DRange *t = dynamic_cast<Param3DRange *>(p);
            Independent3D indep = t->indpendentVariables();

            vec3 val = t->value();

            if (indep == I3D_X)
            {
                all_params[idx] = val.x;
                idx++;
            }
            if (indep == I3D_Y)
            {
                all_params[idx] = val.y;
                idx++;
            }
            if (indep == I3D_Z)
            {
                all_params[idx] = val.z;
                idx++;
            }
            if (indep == I3D_XY)
            {
                all_params[idx]     = val.x;
                all_params[idx + 1] = val.y;
                idx += 2;
            }
            if (indep == I3D_YZ)
            {
                all_params[idx]     = val.y;
                all_params[idx + 1] = val.z;
                idx += 2;
            }
            if (indep == I3D_XZ)
            {
                all_params[idx]     = val.x;
                all_params[idx + 1] = val.z;
                idx += 2;
            }
            if (indep == I3D_XYZ)
            {
                all_params[idx]     = val.x;
                all_params[idx + 1] = val.y;
                all_params[idx + 2] = val.z;
                idx += 3;
            }
        }

        if (p->type == TYPE_RANGE_2D)
        {
            Param2DRange *t = dynamic_cast<Param2DRange *>(p);
            Independent2D indep = t->indpendentVariables();

            vec2 val = t->value();

            if (indep == I2D_X)
            {
                all_params[idx] = val.x;
                idx++;
            }
            if (indep == I2D_Y)
            {
                all_params[idx] = val.y;
                idx++;
            }
            if (indep == I2D_XY)
            {
                all_params[idx]     = val.x;
                all_params[idx + 1] = val.y;
                idx += 2;
            }
        }

        if (p->type == TYPE_POSITION_2D)
        {
            Param2DPosition *t = dynamic_cast<Param2DPosition *>(p);
            Independent2D indep = t->indpendentVariables();

            vec2 val = t->value();

            if (indep == I2D_X)
            {
                all_params[idx] = val.x;
                idx++;
            }
            if (indep == I2D_Y)
            {
                all_params[idx] = val.y;
                idx++;
            }
            if (indep == I2D_XY)
            {
                all_params[idx]     = val.x;
                all_params[idx + 1] = val.y;
                idx += 2;
            }
        }
    }

    return all_params;
}

void Grammar::printParameters()
{
	for (auto iter = m_parameters.begin(); iter != m_parameters.end(); ++iter)
	{
		Param *p = iter.value();		

		if (p->type == TYPE_RANGE_3D)
		{
			Param3DRange *t = dynamic_cast<Param3DRange *>(p);
			//qDebug() << iter.key() << iter.value()->name << t->value().x << t->value().y << t->value().z << t->indpendentVariables();
		}

		if (p->type == TYPE_RANGE_2D)
		{
			Param2DRange *t = dynamic_cast<Param2DRange *>(p);
			//qDebug() << iter.key() << iter.value()->name << t->value().x << t->value().y << t->indpendentVariables();
		}

		if (p->type == TYPE_POSITION_2D)
		{
			Param2DPosition *t = dynamic_cast<Param2DPosition *>(p);
			//qDebug() << iter.key() << iter.value()->name << t->value().x << t->value().y << t->indpendentVariables();
		}
	}
}

int Grammar::determineNrIndependentVariables()
{
	int idx = 0;

    qDebug() << "Size: " << m_parameters.size();
  
	for (auto iter = m_parameters.begin(); iter != m_parameters.end(); ++iter)
	{
		Param *p = iter.value();

        qDebug() << "Parameter ID: " << p->globalId;

		if (p->type == TYPE_RANGE_3D)
		{
			Param3DRange *t = dynamic_cast<Param3DRange *>(p);
			Independent3D indep = t->indpendentVariables();

			if (indep == I3D_X)
			{
				idx++;
			}
			if (indep == I3D_Y)
			{
				idx++;
			}
			if (indep == I3D_Z)
			{
				idx++;
			}
			if (indep == I3D_XY)
			{
				idx += 2;
			}
			if (indep == I3D_YZ)
			{
				idx += 2;
			}
			if (indep == I3D_XZ)
			{
				idx += 2;
			}
			if (indep == I3D_XYZ)
			{
				idx += 3;
			}
		}

		if (p->type == TYPE_RANGE_2D)
		{
			Param2DRange *t = dynamic_cast<Param2DRange *>(p);
			Independent2D indep = t->indpendentVariables();

            qDebug() << indep;

			if (indep == I2D_X)
			{
				idx++;
			}
			if (indep == I2D_Y)
			{
				idx++;
			}
			if (indep == I2D_XY)
			{
				idx += 2;
			}
		}

		if (p->type == TYPE_POSITION_2D)
		{
			Param2DPosition *t = dynamic_cast<Param2DPosition *>(p);
			Independent2D indep = t->indpendentVariables();

			if (indep == I2D_X)
			{
				idx++;
			}
			if (indep == I2D_Y)
			{
				idx++;
			}
			if (indep == I2D_XY)
			{
				idx += 2;
			}
		}
	}

	m_nrIndependentVariables = idx;

	return idx;
}

QMultiMap<QString, Node*> Grammar::nodes()
{
	return m_nodes;
}
