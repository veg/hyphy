#pragma once

#include <QObject>
#include <QNetworkReply>

class QNetworkAccessManager;
class SimpleRequest:public QObject
{
Q_OBJECT
public:
    SimpleRequest();
    ~SimpleRequest(){};
    QNetworkReply* requestUrl(const QString url);
private:
    QNetworkAccessManager* manager;
signals:
    void networkError(const QString);
    void requestFinished();

protected slots:
    void finishReply();
    void errorReply(QNetworkReply::NetworkError);
    void replyFinished(QNetworkReply*);
    const QString parseError(QNetworkReply::NetworkError error);
};
