#include <QNetworkAccessManager>
#include <QDebug>
#include <QFile>
#include <QBuffer>
#include "simplerequest.h"

SimpleRequest::SimpleRequest()
{
        manager = new QNetworkAccessManager;
}

QNetworkReply* SimpleRequest::requestUrl(const QString url)
{
    QNetworkReply* reply = manager->get(QNetworkRequest(QUrl(url)));
    connect(reply, SIGNAL(finished()), this, SLOT(finishReply()));
    return reply;
}

void SimpleRequest::replyFinished(QNetworkReply* reply)
{
    qDebug() << "Request Finished!" << reply;
}

void SimpleRequest::finishReply()
{
    QObject* obj = sender();
    QNetworkReply* reply = qobject_cast<QNetworkReply*>(obj);

    if(reply->error())
    {
       qDebug() << "eror?"; 
        emit networkError(parseError(reply->error()));
        return;
    }

    //emit requestFinished();
}

void SimpleRequest::errorReply(QNetworkReply::NetworkError error)
{
    QString errstr = parseError(error);
    qDebug() << "Error!" << error;
    emit networkError(errstr);
}

const QString SimpleRequest::parseError(QNetworkReply::NetworkError error)
{
    QString errstr = tr("No Error!");
    switch(error)
    {
        case QNetworkReply::ConnectionRefusedError:
        errstr = tr("Connection refused!");
        break;
        case QNetworkReply::HostNotFoundError:
        errstr = tr("Host not found!");
        break;
        case QNetworkReply::RemoteHostClosedError:
        errstr = tr("Remote host closed!");
        break;
        case QNetworkReply::TimeoutError:
        errstr = tr("Timeout!");
        break;
        case QNetworkReply::ContentAccessDenied:
        errstr = tr("Content access denied!");
        break;
        case QNetworkReply::ProtocolFailure:
        errstr = tr("Protocol failure!");
        break;
        case QNetworkReply::ContentNotFoundError:
        errstr = tr("Content not found!");
        break;
        default:
        break;
    }
    //qDebug() << "error" << errstr;
    return errstr;
}
